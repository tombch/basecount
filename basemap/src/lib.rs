use noodles::bam::bai;
// use noodles::bed; // TODO
use noodles::core::region::{Interval, ParseError};
use noodles::core::{Position, Region};
use noodles::sam::alignment::Record;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::sequence::{Base, Sequence};
use noodles::sam::record::{Flags, QualityScores};
use pyo3::exceptions::{PyException, PyIOError};
use pyo3::prelude::*;
use std::cmp;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io;

mod error;
use error::BaseMapError;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy)]
struct Coordinate(usize, usize);

impl IntoPy<PyObject> for Coordinate {
    fn into_py(self, py: Python<'_>) -> PyObject {
        (self.0, self.1).into_py(py)
    }
}

type CoordinateMap = HashMap<Coordinate, [u32; 6]>;

type BaseCount = HashMap<String, CoordinateMap>;

type RefIdToNameMap = HashMap<usize, String>;

type RefNameToLengthMap = HashMap<String, usize>;

const INTEGER_OVERFLOW: &str = "Integer overflow when adding to seq_pos";

fn get_reader(bam_path: String) -> io::Result<noodles::bam::Reader<noodles::bgzf::Reader<File>>> {
    // Open file
    let file = File::open(bam_path)?;

    // Create a reader from the file
    let mut reader = noodles::bam::Reader::new(file);

    // Read the SAM header
    reader.read_header()?;

    // Return the reader
    Ok(reader)
}

fn count_base(
    ref_map: &mut CoordinateMap,
    seq: &Sequence,
    ref_pos: usize,
    seq_pos: Position,
    ins_pos: usize,
) {
    // Match the base at the given seq_pos, and update the CoordinateMap
    match seq.get(seq_pos) {
        Some(&Base::A) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert([0; 6])[0] += 1
        }
        Some(&Base::C) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert([0; 6])[1] += 1
        }
        Some(&Base::G) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert([0; 6])[2] += 1
        }
        Some(&Base::T) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert([0; 6])[3] += 1
        }
        Some(&Base::N) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert([0; 6])[5] += 1
        }
        Some(x) => panic!("Encountered invalid read base: {}", x),
        None => panic!("Attempted to access invalid read index: {}", seq_pos.get()),
    }
}

fn count_record(
    ref_map: &mut CoordinateMap,
    record: &Record,
    base_quality: usize,
) -> Result<(), Box<dyn Error>> {
    // Positions are 1-based
    // This is the start position of the read in the reference
    let mut ref_pos = record
        .alignment_start()
        .ok_or("Could not find alignment_start for record")?
        .get();

    // This is the position locally along the sequence (minimum is 1)
    let mut seq_pos = Position::MIN;

    // The read sequence
    let seq = record.sequence();

    // The read sequence quality scores
    let quals = record.quality_scores();

    // Iterate through CIGAR information
    for cig in record.cigar().iter() {
        match cig.kind() {
            // Match/mismatch consumes both the reference and sequence
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 1..=cig.len() {
                    if min_base_quality(quals, seq_pos, base_quality)? {
                        count_base(ref_map, seq, ref_pos, seq_pos, 0);
                    }

                    ref_pos += 1;
                    seq_pos = seq_pos.checked_add(1).ok_or(INTEGER_OVERFLOW)?;
                }
            }

            // Insertion consumes the sequence only
            Kind::Insertion => {
                for i in 1..=cig.len() {
                    if min_base_quality(quals, seq_pos, base_quality)? {
                        count_base(ref_map, seq, ref_pos, seq_pos, i);
                    }

                    seq_pos = seq_pos.checked_add(1).ok_or(INTEGER_OVERFLOW)?;
                }
            }

            // Deletion/skip consumes the reference only
            Kind::Deletion | Kind::Skip => {
                for _ in 1..=cig.len() {
                    ref_map.entry(Coordinate(ref_pos, 0)).or_insert([0; 6])[4] += 1;
                    ref_pos += 1;
                }
            }

            // Softclip consumes the sequence only
            Kind::SoftClip => {
                seq_pos = seq_pos.checked_add(cig.len()).ok_or(INTEGER_OVERFLOW)?;
            }

            // Hardclip and padding doesn't consume the reference or sequence
            Kind::HardClip | Kind::Pad => {}
        };
    }
    Ok(())
}

fn init_maps() -> (RefIdToNameMap, RefNameToLengthMap, BaseCount) {
    // Map of reference sequence ids to their reference name
    let ref_seq_ids: RefIdToNameMap = HashMap::new();

    // Map of reference names to reference lengths
    let ref_lengths: RefNameToLengthMap = HashMap::new();

    // Map of maps, storing count data for each reference
    let map: BaseCount = BaseCount::new();

    (ref_seq_ids, ref_lengths, map)
}

fn init_region_coordinates(
    map: &mut BaseCount,
    region: &Region,
    ref_lengths: HashMap<String, usize>,
) -> Result<(), Box<dyn Error>> {
    let region_name = region.name();
    let interval = region.interval();

    // Get length of the region name's sequence
    let ref_length = ref_lengths
        .get(region_name)
        .ok_or("Could not find region name in references")?
        .to_owned();

    // Handle unbounded region start
    let region_start = match interval.start() {
        Some(x) => x.get(),
        None => 1,
    };

    // Handle unbounded region end
    let region_end = match interval.end() {
        Some(x) => x.get(),
        None => ref_length,
    };

    // Insert Coordinates for all region positions
    for i in cmp::min(region_start, ref_length)..=cmp::min(region_end, ref_length) {
        map.entry(region_name.to_owned())
            .or_default()
            .entry(Coordinate(i, 0))
            .or_insert([0; 6]);
    }

    Ok(())
}

fn record_in_region(
    record: &Record,
    ref_seq_ids: &HashMap<usize, String>,
    region: &Region,
) -> Result<bool, Box<dyn Error>> {
    let record_ref = ref_seq_ids
        .get(
            &record
                .reference_sequence_id()
                .expect("Could not read reference_sequence_id"),
        )
        .expect("Could not get reference_seq_name using reference_seq_id");

    let seq_start = record
        .alignment_start()
        .ok_or("Could not read sequence start position")?;

    let seq_end = record
        .alignment_end()
        .ok_or("Could not read sequence end position")?;

    let seq_interval = Interval::from(seq_start..=seq_end);

    if record_ref == region.name() && region.interval().intersects(seq_interval) {
        Ok(true)
    } else {
        Ok(false)
    }
}

fn min_base_quality(
    quals: &QualityScores,
    seq_pos: Position,
    base_quality: usize,
) -> Result<bool, Box<dyn Error>> {
    let base_qual = usize::from(
        quals
            .get(seq_pos)
            .ok_or("Could not get quality_score for base")?
            .get(),
    );

    if base_qual >= base_quality {
        Ok(true)
    } else {
        Ok(false)
    }
}

fn min_mapping_quality(record: &Record, mapping_quality: usize) -> Result<bool, Box<dyn Error>> {
    let map_qual = usize::from(
        record
            .mapping_quality()
            .ok_or("Could not read mapping_quality for record")?
            .get(),
    );

    if map_qual >= mapping_quality {
        Ok(true)
    } else {
        Ok(false)
    }
}

#[pyfunction]
fn all_(bam_path: String, mapping_quality: usize, base_quality: usize) -> PyResult<BaseCount> {
    // Map of maps, storing count data for each reference
    let mut map: BaseCount = BaseCount::new();

    // Reader for iterating through records
    let mut reader = get_reader(bam_path).map_err(|x| PyIOError::new_err(x.to_string()))?;

    // Reference sequence information
    let ref_seqs = reader.read_reference_sequences()?;

    // Add each reference to the map
    // For each reference, insert Coordinates for all reference positions
    for reff in ref_seqs.iter() {
        let ref_name = reff.0.to_owned();
        let ref_length = reff.1.length().get();

        for i in 1..=ref_length {
            map.entry(ref_name.to_owned())
                .or_default()
                .entry(Coordinate(i, 0))
                .or_insert([0; 6]);
        }
    }

    let flags = Flags::from(
        Flags::UNMAPPED.bits()
            + Flags::SUPPLEMENTARY.bits()
            + Flags::SECONDARY.bits()
            + Flags::QC_FAIL.bits()
            + Flags::DUPLICATE.bits(),
    );

    for result in reader.records() {
        let record = result?;

        if record.flags().intersects(flags)
            || !min_mapping_quality(&record, mapping_quality)
                .map_err(|x| PyException::new_err(x.to_string()))?
        {
            continue;
        }

        let ref_seq_id = record
            .reference_sequence_id()
            .expect("Could not read reference_sequence_id for record");

        let ref_name = ref_seqs
            .get_index(ref_seq_id)
            .expect("Could not get reference_seq_name using reference_sequence_id")
            .0;

        let ref_map = map
            .get_mut(ref_name)
            .expect("Could not find reference_seq_name in map");

        count_record(ref_map, &record, base_quality)
            .map_err(|x| PyException::new_err(x.to_string()))?;
    }
    Ok(map)
}

#[pyfunction]
fn query_(
    bam_path: String,
    region: String,
    mapping_quality: usize,
    base_quality: usize,
) -> PyResult<BaseCount> {
    // Create initial maps
    let (mut ref_seq_ids, mut ref_lengths, mut map) = init_maps();

    // Reader for iterating through records
    let mut reader = get_reader(bam_path).map_err(|x| PyIOError::new_err(x.to_string()))?;

    // Reference sequence information
    let ref_seqs = reader.read_reference_sequences()?;

    // Add reference sequence information to HashMaps
    for (i, reff) in ref_seqs.iter().enumerate() {
        ref_seq_ids.insert(i, reff.0.to_owned());
        ref_lengths.insert(reff.0.to_owned(), reff.1.length().get());
    }

    // Parse region
    let region: Region = region
        .parse()
        .map_err(|x: ParseError| PyException::new_err(x.to_string()))?;
    let region_name = region.name();

    // Create coordinates for all points defined in the region
    init_region_coordinates(&mut map, &region, ref_lengths)
        .map_err(|x| PyException::new_err(x.to_string()))?;

    let ref_map = map
        .get_mut(region_name)
        .expect("Could not find region_name in map");

    let flags = Flags::from(
        Flags::UNMAPPED.bits()
            + Flags::SUPPLEMENTARY.bits()
            + Flags::SECONDARY.bits()
            + Flags::QC_FAIL.bits()
            + Flags::DUPLICATE.bits(),
    );

    for result in reader.records() {
        let record = result?;

        if record.flags().intersects(flags)
            || !record_in_region(&record, &ref_seq_ids, &region)
                .map_err(|x| PyException::new_err(x.to_string()))?
            || !min_mapping_quality(&record, mapping_quality)
                .map_err(|x| PyException::new_err(x.to_string()))?
        {
            continue;
        }

        count_record(ref_map, &record, base_quality)
            .map_err(|x| PyException::new_err(x.to_string()))?;
    }
    Ok(map)
}

#[pyfunction]
fn iquery_(
    bam_path: String,
    bai_path: String,
    region: String,
    mapping_quality: usize,
    base_quality: usize,
) -> PyResult<BaseCount> {
    // Create initial maps
    let (mut ref_seq_ids, mut ref_lengths, mut map) = init_maps();

    // Reader for iterating through records
    let mut reader = get_reader(bam_path).map_err(|x| PyIOError::new_err(x.to_string()))?;

    // Reference sequence information
    let ref_seqs = reader
        .read_reference_sequences()
        .map_err(|x| PyException::new_err(x.to_string()))?;

    // Add reference sequence information to HashMaps
    for (i, reff) in ref_seqs.iter().enumerate() {
        ref_seq_ids.insert(i, reff.0.to_owned());
        ref_lengths.insert(reff.0.to_owned(), reff.1.length().get());
    }

    // Parse region
    let region: Region = region
        .parse()
        .map_err(|x: ParseError| PyException::new_err(x.to_string()))?;
    let region_name = region.name();

    // Read the index file
    let index = bai::read(bai_path)?;

    // Create query iterator over reads intersecting the region
    let query = reader.query(&ref_seqs, &index, &region)?;

    // Create coordinates for all points defined in the region
    init_region_coordinates(&mut map, &region, ref_lengths)
        .map_err(|x| PyException::new_err(x.to_string()))?;

    let ref_map = map
        .get_mut(region_name)
        .expect("Could not find region_name in map");

    let flags = Flags::from(
        Flags::UNMAPPED.bits()
            + Flags::SUPPLEMENTARY.bits()
            + Flags::SECONDARY.bits()
            + Flags::QC_FAIL.bits()
            + Flags::DUPLICATE.bits(),
    );

    for result in query {
        let record = result?;

        if record.flags().intersects(flags)
            || !record_in_region(&record, &ref_seq_ids, &region)
                .map_err(|x| PyException::new_err(x.to_string()))?
            || !min_mapping_quality(&record, mapping_quality)
                .map_err(|x| PyException::new_err(x.to_string()))?
        {
            continue;
        }

        count_record(ref_map, &record, base_quality)
            .map_err(|x| PyException::new_err(x.to_string()))?;
    }
    Ok(map)
}

#[pyfunction]
fn parse_region_(region: String) -> PyResult<(String, usize, usize)> {
    let region: Region = region
        .parse()
        .map_err(|x: ParseError| PyException::new_err(x.to_string()))?;
    let interval = region.interval();
    let start = match interval.start() {
        Some(x) => x.get(),
        None => 0,
    };
    let end = match interval.end() {
        Some(x) => x.get(),
        None => 0,
    };

    Ok((region.name().to_string(), start, end))
}

/// A Python module implemented in Rust.
#[pymodule]
fn basemap(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(all_, m)?)?;
    m.add_function(wrap_pyfunction!(query_, m)?)?;
    m.add_function(wrap_pyfunction!(iquery_, m)?)?;
    m.add_function(wrap_pyfunction!(parse_region_, m)?)?;

    Ok(())
}
