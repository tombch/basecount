use noodles::bam::bai;
use noodles::core::region::{Interval, ParseError};
use noodles::core::{Position, Region};
use noodles::sam::alignment::Record;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::sequence::{Base, Sequence};
use noodles::sam::record::{Flags, QualityScores};
use pyo3::exceptions::{PyException, PyIOError, PyIndexError, PyKeyError, PyOverflowError};
use pyo3::prelude::*;
use std::collections::HashMap;
use std::fs::File;

// use noodles::bed; // TODO

mod error;
use error::BaseMapError;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy)]
struct Coordinate(usize, usize);

impl IntoPy<PyObject> for Coordinate {
    fn into_py(self, py: Python<'_>) -> PyObject {
        (self.0, self.1).into_py(py)
    }
}

impl From<BaseMapError> for PyErr {
    fn from(e: BaseMapError) -> Self {
        match e {
            BaseMapError::KeyNotFound => PyKeyError::new_err(e.to_string()),
            BaseMapError::IndexNotFound => PyIndexError::new_err(e.to_string()),
            BaseMapError::IntegerOverflow => PyOverflowError::new_err(e.to_string()),
            BaseMapError::IOError(e) => PyIOError::new_err(e.to_string()),
            _ => PyException::new_err(e.to_string()),
        }
    }
}

type CoordinateMap = HashMap<Coordinate, [u32; 6]>;

type BaseMap = HashMap<String, CoordinateMap>;

type RefLengths = HashMap<String, usize>;

/// Open the BAM file located at `bam_path` and return a reader.
fn get_reader(
    bam_path: String,
) -> Result<noodles::bam::Reader<noodles::bgzf::Reader<File>>, BaseMapError> {
    // Open file
    let file = File::open(bam_path)?;

    // Create a reader from the file
    let mut reader = noodles::bam::Reader::new(file);

    // Read the SAM header
    reader.read_header()?;

    // Return the reader
    Ok(reader)
}

/// Add the base from `seq` at `(seq_pos, ins_pos)` to `ref_map`.
fn count_base(
    ref_map: &mut CoordinateMap,
    seq: &Sequence,
    ref_pos: usize,
    seq_pos: Position,
    ins_pos: usize,
) -> Result<(), BaseMapError> {
    // Match the base at the given seq_pos, and update the CoordinateMap
    match seq.get(seq_pos) {
        Some(&Base::A) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[0] += 1;
            Ok(())
        }
        Some(&Base::C) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[1] += 1;
            Ok(())
        }
        Some(&Base::G) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[2] += 1;
            Ok(())
        }
        Some(&Base::T) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[3] += 1;
            Ok(())
        }
        Some(&Base::N) => {
            ref_map
                .entry(Coordinate(ref_pos, ins_pos))
                .or_insert_with(|| [0; 6])[5] += 1;
            Ok(())
        }
        Some(_) => Err(BaseMapError::InvalidBase),
        None => Err(BaseMapError::KeyNotFound),
    }
}

/// Use the CIGAR information of `record` to count each base in its sequence, and add them to `ref_map`.
///
/// Bases are ignored if their quality score is less than `base_quality`.
fn count_record(
    ref_map: &mut CoordinateMap,
    record: &Record,
    base_quality: usize,
) -> Result<(), BaseMapError> {
    // Positions are 1-based
    // This is the start position of the read in the reference
    let mut ref_pos = record
        .alignment_start()
        .ok_or_else(|| BaseMapError::AlignmentStartNotFound)?
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
                        count_base(ref_map, seq, ref_pos, seq_pos, 0)?;
                    }

                    ref_pos += 1;
                    seq_pos = seq_pos
                        .checked_add(1)
                        .ok_or_else(|| BaseMapError::IntegerOverflow)?;
                }
            }

            // Insertion consumes the sequence only
            Kind::Insertion => {
                for i in 1..=cig.len() {
                    if min_base_quality(quals, seq_pos, base_quality)? {
                        count_base(ref_map, seq, ref_pos, seq_pos, i)?;
                    }

                    seq_pos = seq_pos
                        .checked_add(1)
                        .ok_or_else(|| BaseMapError::IntegerOverflow)?;
                }
            }

            // Deletion/skip consumes the reference only
            Kind::Deletion | Kind::Skip => {
                for _ in 1..=cig.len() {
                    ref_map
                        .entry(Coordinate(ref_pos, 0))
                        .or_insert_with(|| [0; 6])[4] += 1;

                    ref_pos += 1;
                }
            }

            // Softclip consumes the sequence only
            Kind::SoftClip => {
                seq_pos = seq_pos
                    .checked_add(cig.len())
                    .ok_or_else(|| BaseMapError::IntegerOverflow)?;
            }

            // Hardclip and padding don't consume the reference or the sequence
            Kind::HardClip | Kind::Pad => {}
        };
    }
    Ok(())
}

fn init_maps() -> (RefLengths, BaseMap) {
    // Map of reference names to reference lengths
    let ref_lengths: RefLengths = HashMap::new();

    // Map of maps, storing count data for each reference
    let map: BaseMap = BaseMap::new();

    (ref_lengths, map)
}

fn init_region_coordinates(
    map: &mut BaseMap,
    region: &Region,
    ref_lengths: HashMap<String, usize>,
) -> Result<(), BaseMapError> {
    let region_name = region.name();
    let interval = region.interval();

    // Get length of the region name's sequence
    let ref_length = ref_lengths
        .get(region_name)
        .ok_or_else(|| BaseMapError::KeyNotFound)?
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

    // Create coordinate map for the particular region
    let ref_map = map.entry(region_name.to_owned()).or_default();

    // Insert Coordinates for all region positions
    for i in region_start..=region_end {
        ref_map.entry(Coordinate(i, 0)).or_insert_with(|| [0; 6]);
    }

    Ok(())
}

/// Check the interval defined by the alignment of `record` intersects the interval defined in `region`.
fn intersects(record: &Record, region: &Region) -> Result<bool, BaseMapError> {
    let seq_start = record
        .alignment_start()
        .ok_or_else(|| BaseMapError::AlignmentStartNotFound)?;

    let seq_end = record
        .alignment_end()
        .ok_or_else(|| BaseMapError::AlignmentEndNotFound)?;

    let seq_interval = Interval::from(seq_start..=seq_end);

    if region.interval().intersects(seq_interval) {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Check the quality score for the base at `seq_pos` is greater than or equal to `base_quality`.
fn min_base_quality(
    quals: &QualityScores,
    seq_pos: Position,
    base_quality: usize,
) -> Result<bool, BaseMapError> {
    let base_qual = usize::from(
        quals
            .get(seq_pos)
            .ok_or_else(|| BaseMapError::QualityScoreNotFound)?
            .get(),
    );

    if base_qual >= base_quality {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Check the mapping score for `record` is greater than or equal to `mapping_quality`.
fn min_mapping_quality(record: &Record, mapping_quality: usize) -> Result<bool, BaseMapError> {
    let map_qual = usize::from(
        record
            .mapping_quality()
            .ok_or_else(|| BaseMapError::MappingQualityNotFound)?
            .get(),
    );

    if map_qual >= mapping_quality {
        Ok(true)
    } else {
        Ok(false)
    }
}

#[pyfunction]
fn all_(bam_path: String, mapping_quality: usize, base_quality: usize) -> PyResult<BaseMap> {
    // Create initial maps
    let (mut ref_lengths, mut map) = init_maps();

    // Reader for iterating through records
    let mut reader = get_reader(bam_path)?;

    // Reference sequence information
    let ref_seqs = reader.read_reference_sequences()?;

    // Add reference sequence information to HashMaps
    for reff in ref_seqs.iter() {
        ref_lengths.insert(reff.0.to_owned(), reff.1.length().get());
    }

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

        if record.flags().intersects(flags) || !min_mapping_quality(&record, mapping_quality)? {
            continue;
        }

        let ref_seq_id = record
            .reference_sequence_id()
            .ok_or_else(|| BaseMapError::ReferenceSequenceIDNotFound)?;

        let ref_name = ref_seqs
            .get_index(ref_seq_id)
            .ok_or_else(|| BaseMapError::KeyNotFound)?
            .0;

        let ref_map = map
            .get_mut(ref_name)
            .ok_or_else(|| BaseMapError::KeyNotFound)?;

        count_record(ref_map, &record, base_quality)?;
    }
    Ok(map)
}

#[pyfunction]
fn query_(
    bam_path: String,
    bai_path: Option<String>,
    region: String,
    mapping_quality: usize,
    base_quality: usize,
) -> PyResult<BaseMap> {
    // Create initial maps
    let (mut ref_lengths, mut map) = init_maps();

    // Reader for iterating through records
    let mut reader = get_reader(bam_path)?;

    // Reference sequence information
    let ref_seqs = reader.read_reference_sequences()?;

    // Add reference sequence information to HashMaps
    for reff in ref_seqs.iter() {
        ref_lengths.insert(reff.0.to_owned(), reff.1.length().get());
    }

    // Parse region
    let region: Region = region
        .parse()
        .map_err(|x: ParseError| PyException::new_err(x.to_string()))?;
    let region_name = region.name();

    // Create coordinates for all points defined in the region
    init_region_coordinates(&mut map, &region, ref_lengths)?;

    let flags = Flags::from(
        Flags::UNMAPPED.bits()
            + Flags::SUPPLEMENTARY.bits()
            + Flags::SECONDARY.bits()
            + Flags::QC_FAIL.bits()
            + Flags::DUPLICATE.bits(),
    );

    // Get Coordinate map for particular region
    let ref_map = map
        .get_mut(region_name)
        .ok_or_else(|| BaseMapError::KeyNotFound)?;

    if let Some(b_path) = bai_path {
        // Read the index file
        let index = bai::read(b_path)?;

        // Create query iterator over reads intersecting the region
        let query = reader.query(&ref_seqs, &index, &region)?;

        for result in query {
            let record = result?;
            if record.flags().intersects(flags) || !min_mapping_quality(&record, mapping_quality)? {
                continue;
            }

            count_record(ref_map, &record, base_quality)?;
        }
    } else {
        for result in reader.records() {
            let record = result?;
            let record_ref_name = ref_seqs
                .get_index(
                    record
                        .reference_sequence_id()
                        .ok_or_else(|| BaseMapError::ReferenceSequenceIDNotFound)?,
                )
                .ok_or_else(|| BaseMapError::IndexNotFound)?
                .0;

            if record.flags().intersects(flags)
                || record_ref_name != region.name()
                || !intersects(&record, &region)?
                || !min_mapping_quality(&record, mapping_quality)?
            {
                continue;
            }

            count_record(ref_map, &record, base_quality)?;
        }
    }

    Ok(map)
}

#[pyfunction]
fn parse_region_(region: String) -> PyResult<(String, Option<usize>, Option<usize>)> {
    let region: Region = region
        .parse()
        .map_err(|x: ParseError| PyException::new_err(x.to_string()))?;
    let interval = region.interval();
    let start = match interval.start() {
        Some(x) => Some(x.get()),
        None => None,
    };
    let end = match interval.end() {
        Some(x) => Some(x.get()),
        None => None,
    };

    Ok((region.name().to_string(), start, end))
}

/// A Python module implemented in Rust.
#[pymodule]
fn basemap(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(all_, m)?)?;
    m.add_function(wrap_pyfunction!(query_, m)?)?;
    m.add_function(wrap_pyfunction!(parse_region_, m)?)?;

    Ok(())
}
