import os
import math
import glob
import pysam
import pytest
import itertools
import numpy as np
import concurrent.futures
from count import bcount
from basecount import BaseCount


# Put path to directory of BAM files here
bams_dir = "/home/tom/git/samples/"

# List of files within the directory that end in .bam
bams = glob.glob(f"{bams_dir}/*.bam")

# Min base quality and min mapping quality test parameters
min_base_qualities = [0, 10, 20, 30, 40]
min_mapping_qualities = [0, 10, 20, 30, 40, 50, 60]

# Test params assembled together
params = list(itertools.product(bams, min_base_qualities, min_mapping_qualities))

# Max number of workers for getting test data
num_workers = 8


def calculate_ref_region(bam, ref, start, end, min_base_quality, min_mapping_quality):
    samfile = pysam.AlignmentFile(bam, mode="rb")
    normalising_factor = (1 / math.log2(6)) # Normalises maximum entropy to 1
    secondary_normalising_factor = (1 / math.log2(5)) # Normalises maximum secondary_entropy to 1
    empty_bases = [0, 0, 0, 0, 0, 0]
    invalid_base_percentages = [-1, -1, -1, -1, -1, -1] # Percentage values where coverage is zero
    
    base_count = [{'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'DS' : 0, 'N' : 0} for _ in range(end - start)]
    region = [[] for _ in range(end - start)]

    for pileupcolumn in samfile.pileup(ref, start, end, min_base_quality=0, min_mapping_quality=0, max_depth=100000000, stepper="nofilter"):
        ref_pos = pileupcolumn.reference_pos # type: ignore
        if start <= ref_pos < end:
            for pileupread in pileupcolumn.pileups: # type: ignore
                if (pileupread.alignment.mapping_quality < min_mapping_quality):
                    continue    
                if not pileupread.is_del and not pileupread.is_refskip:
                    if pileupread.alignment.query_qualities[pileupread.query_position] >= min_base_quality:
                        base_count[ref_pos - start][pileupread.alignment.query_sequence[pileupread.query_position]] += 1
                else:
                    base_count[ref_pos - start]['DS'] += 1

            # TODO: Had to change this to fit with the min base and mapping args, which begs question how should coverage be counted ?
            # coverage = pileupcolumn.nsegments # type: ignore
            coverage = sum(base_count[ref_pos - start].values())

            if coverage != 0:
                b_c = list(base_count[ref_pos - start].values())

                base_probabilities = [count / coverage for count in b_c]
                base_percentages = [100 * probability for probability in base_probabilities]
                entropy = normalising_factor * sum([-(x * math.log2(x)) if x != 0 else 0 for x in base_probabilities])

                secondary_base_count = list(b_c)
                secondary_base_count.pop(np.argmax(b_c))
                secondary_coverage = sum(secondary_base_count)
                if secondary_coverage != 0:
                    secondary_base_probabilities = [count / secondary_coverage for count in secondary_base_count]
                    secondary_entropy = secondary_normalising_factor * sum([-(x * math.log2(x)) if x != 0 else 0 for x in secondary_base_probabilities])
                else:
                    secondary_entropy = 1
            else:
                base_percentages = invalid_base_percentages
                entropy = 1
                secondary_entropy = 1

            region[ref_pos - start].append(ref)
            region[ref_pos - start].append(ref_pos + 1)
            region[ref_pos - start].append(coverage)
            region[ref_pos - start].extend(base_count[ref_pos - start].values())
            region[ref_pos - start].extend(base_percentages)
            region[ref_pos - start].append(entropy)
            region[ref_pos - start].append(secondary_entropy)
    
    # Fill empty rows
    for i, row in enumerate(region):
        if len(row) == 0:
            row.append(ref)
            row.append(start + i + 1)
            row.append(0)
            row.extend(empty_bases)
            row.extend(invalid_base_percentages)
            row.extend([1, 1])

    samfile.close()
    return region


def get_test_data(bam, min_base_quality=0, min_mapping_quality=0):    
    samfile = pysam.AlignmentFile(bam, "rb")
    data = {}

    for ref in samfile.references:
        ref_len = samfile.lengths[samfile.references.index(ref)]
        num_reads = len([read for read in samfile.fetch(ref) if read.mapping_quality >= min_mapping_quality])
        region_markers = np.linspace(0, ref_len, num=num_workers + 1, dtype=int)
        regions = [(region_markers[i], region_markers[i + 1]) for i in range(len(region_markers) - 1)]

        # Count the bases
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = {start : executor.submit(calculate_ref_region, bam, ref, start, end, min_base_quality, min_mapping_quality) for start, end in regions}

        ref_data = []
        for start, _ in regions:
            ref_data.extend(results[start].result())
        data[ref] = ref_data, num_reads

    samfile.close()
    return data


def get_entropy(probabilities):
    return sum([-(x * math.log2(x)) if x != 0 else 0 for x in probabilities])


def get_stats(base_counts, ref, show_n_bases=False, long_format=False):
    # Order of bases in each output row of the basecount
    bases = ["A", "C", "G", "T", "DS", "N"]
    n_index = bases.index("N")

    if not show_n_bases:
        bases.pop(n_index)

    num_bases = len(bases)
    invalid_percentages = [-1] * num_bases
    normalising_factor = 1 / math.log2(num_bases)
    secondary_normalising_factor = 1 / math.log2(num_bases - 1)

    # Generate per-position statistics
    data = []
    for reference_pos, base_count in enumerate(base_counts):
        if not show_n_bases:
            base_count.pop(n_index)
        
        # Defaults for when the coverage is zero
        base_percentages = invalid_percentages
        entropy = 1
        secondary_entropy = 1
        coverage = sum(base_count)

        if coverage != 0:
            base_probabilities = [count / coverage for count in base_count]
            base_percentages = [100 * probability for probability in base_probabilities]
            entropy = normalising_factor * get_entropy(base_probabilities)

            secondary_base_count = list(base_count)
            secondary_base_count.pop(np.argmax(base_count))
            secondary_coverage = sum(secondary_base_count)
            if secondary_coverage != 0:
                secondary_base_probabilities = [count / secondary_coverage for count in secondary_base_count]
                secondary_entropy = secondary_normalising_factor * get_entropy(secondary_base_probabilities)

        # Create row of basecount data, either in long format or wide format 
        # Wide format is default
        if long_format:
            for base, count, percentage in zip(bases, base_count, base_percentages):
                row = []
                row.append(ref)
                row.append(reference_pos + 1) # Output one-based coordinates
                row.append(coverage)
                row.append(base)
                row.append(count)
                row.append(percentage)
                row.append(entropy)
                row.append(secondary_entropy)
                data.append(row)
        else:
            row = []
            row.append(ref)
            row.append(reference_pos + 1) # Output one-based coordinates
            row.append(coverage)
            row.extend(base_count)
            row.extend(base_percentages)
            row.append(entropy)
            row.append(secondary_entropy)
            data.append(row)
    
    return data


def get_references(samfile, references=None):
    # Determine which references will be covered
    if references is None:
        # Calculate data on all references by default
        references = samfile.references
    else:
        # If references are provided, determine their validity
        for reference in references:
            if not (reference in samfile.references):
                raise Exception(f"{reference} is not a valid reference")
    return set(references)


def open_samfile(bam):
    # Temporarily suppress HTSlib's messages when opening the file
    old_verbosity = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(bam, mode="rb")
    pysam.set_verbosity(old_verbosity)
    return samfile


def get_basecounts(bam, references=None, min_base_quality=0, min_mapping_quality=0, show_n_bases=False, long_format=False):
    samfile = open_samfile(bam)
    references = get_references(samfile, references)

    # Iterator through all reads in the file (including unmapped)
    all_reads = samfile.fetch(until_eof=True) 
    
    # Prepare read_data structure
    read_data = {}
    for ref in references:
        read_data[ref] = {}
        read_data[ref]["ref_len"] = samfile.lengths[samfile.references.index(ref)]
        read_data[ref]["reads"] = []
        read_data[ref]["qualities"] = []
        read_data[ref]["starts"] = []
        read_data[ref]["ctuples"] = []

    # Filter for mapped reads, storing each read's sequence, start and ctuples
    for read in all_reads:
        if (not read.is_unmapped) and (read.mapping_quality >= min_mapping_quality):
            ref = read.reference_name
            if ref in references:
                read_data[ref]["reads"].append(read.query_alignment_sequence)
                read_data[ref]["qualities"].append(read.query_alignment_qualities)
                read_data[ref]["starts"].append(read.reference_start)
                read_data[ref]["ctuples"].append(read.cigartuples)
        
    # For each reference, verify the read data and calculate num_reads
    for ref, ref_data in read_data.items():
        if len({len(ref_data[x]) for x in ["reads", "qualities", "starts", "ctuples"]}) != 1:
            raise Exception(f"The number of reads, qualities, starts and ctuples for {ref} do not match")
        ref_data["num_reads"] = len(ref_data["reads"])

    basecount_data = {}    
    for ref in references:
        # Count the bases
        base_counts = bcount(
            read_data[ref]["ref_len"],
            min_base_quality,
            read_data[ref]["reads"],
            read_data[ref]["qualities"],
            read_data[ref]["starts"],
            read_data[ref]["ctuples"]
        )
            
        basecount_data[ref] = {
            "rows" : get_stats(
                    base_counts,
                    ref,
                    show_n_bases=show_n_bases,
                    long_format=long_format
                ), 
            "num_reads" : read_data[ref]["num_reads"]
        }

    samfile.close()
    return basecount_data


# @pytest.mark.parametrize("bam,mbq,mmq", params)
# def test_basecount(bam, mbq, mmq):
#     # Create an index (if it doesn't already exist) in the same dir as the BAM
#     if not os.path.isfile(bam + '.bai'):
#         pysam.index(bam) # type: ignore
    
#     bc = BaseCount(bam, min_base_quality=mbq, min_mapping_quality=mmq, show_n_bases=True)
#     test_data = get_test_data(bam, min_base_quality=mbq, min_mapping_quality=mmq)

#     # Test for matching references
#     assert bc.references == list(test_data.keys())

#     # Compare the data for each reference
#     for ref in bc.references:
#         test_ref_data, test_ref_num_reads = test_data[ref]

#         # Test for matching reference lengths
#         assert bc.reference_lengths[ref] == len(test_ref_data)

#         # Test for matching total of reads
#         # Given the other tests, and how this value is calculated, this test is a bit pointless
#         # But no reason not to include it for completion
#         assert bc.num_reads(ref) == test_ref_num_reads

#         # Iterate through rows, comparing each
#         for bc_row, test_row in zip(bc.rows(ref), test_ref_data):
#             assert bc_row == test_row


@pytest.mark.parametrize("bam,mbq,mmq", params)
def test_basecount_old(bam, mbq, mmq):
    # Create an index (if it doesn't already exist) in the same dir as the BAM
    if not os.path.isfile(bam + '.bai'):
        pysam.index(bam) # type: ignore
    
    bc = BaseCount(bam, min_base_quality=mbq, min_mapping_quality=mmq, show_n_bases=True)
    bc_old = get_basecounts(bam, references=None, min_base_quality=mbq, min_mapping_quality=mmq, show_n_bases=True)

    # Test for matching references
    assert bc.references == list(bc_old.keys())

    # Compare the data for each reference
    for ref in bc.references:
        # Test for matching reference lengths
        assert bc.reference_lengths[ref] == len(bc_old[ref]["rows"])

        # Test for matching total of reads
        # Given the other tests, and how this value is calculated, this test is a bit pointless
        # But no reason not to include it for completion
        assert bc.num_reads(ref) == bc_old[ref]["num_reads"]

        # Iterate through rows, comparing each
        for bc_row, test_row in zip(bc.rows(ref), bc_old[ref]["rows"]):
            assert bc_row == test_row
