import os
import math
import glob
import pysam
import pytest
import numpy as np
import concurrent.futures
from basecount import BaseCount


# Put path to directory of BAM files here
bams_dir = ""
# List of files within the directory that end in .bam
bams = glob.glob(f'{bams_dir}/*.bam')
# Max number of workers for getting test data
num_workers = 8


def calculate_ref_region(bam, ref, start, end):
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
                if not pileupread.is_del and not pileupread.is_refskip:
                    base_count[ref_pos - start][pileupread.alignment.query_sequence[pileupread.query_position]] += 1
                else:
                    base_count[ref_pos - start]['DS'] += 1

            coverage = pileupcolumn.nsegments # type: ignore
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


def get_test_data(bam):    
    samfile = pysam.AlignmentFile(bam, "rb")
    data = {}

    for ref in samfile.references:
        ref_len = samfile.lengths[samfile.references.index(ref)]
        num_reads = len([read for read in samfile.fetch(ref)])
        region_markers = np.linspace(0, ref_len, num=num_workers + 1, dtype=int)
        regions = [(region_markers[i], region_markers[i + 1]) for i in range(len(region_markers) - 1)]

        # Count the bases
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = {start : executor.submit(calculate_ref_region, bam, ref, start, end) for start, end in regions}

        ref_data = []
        for start, _ in regions:
            ref_data.extend(results[start].result())
        data[ref] = ref_data, num_reads

    samfile.close()
    return data


@pytest.mark.parametrize("bam", bams)
def test_basecount(bam):
    # Create an index (if it doesn't already exist) in the same dir as the BAM
    if not os.path.isfile(bam + '.bai'):
        pysam.index(bam) # type: ignore
    
    bc = BaseCount(bam)
    test_data = get_test_data(bam)

    # Test for matching references
    assert bc.references == list(test_data.keys())

    # Compare the data for each reference
    for ref in bc.references:
        test_ref_data, test_ref_num_reads = test_data[ref]

        # Test for matching reference lengths
        assert bc.reference_lengths[ref] == len(test_ref_data)

        # Test for matching total of reads
        assert bc.num_reads(ref) == test_ref_num_reads

        # Iterate through rows, comparing each
        for bc_row, test_row in zip(bc.rows(ref), test_ref_data):
            assert bc_row == test_row
