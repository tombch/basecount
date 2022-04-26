import os
import math
import glob
import pysam
import pytest
import basecount


# Put path to directory of BAM files here
bams_dir = ""
# List of files within the directory that end in .bam
bams = glob.glob(f'{bams_dir}/*.bam')


def get_test_data(bam):    
    samfile = pysam.AlignmentFile(bam, "rb")

    data = {}
    normalising_factor = 1 / math.log2(5) # Normalises maximum entropy to 1
    invalid_base_percentages = [-1, -1, -1, -1, -1, -1] # Percentage values where coverage is zero

    for ref in samfile.references:
        num_reads = len([read for read in samfile.fetch(ref)])
        ref_len = samfile.lengths[samfile.references.index(ref)]
        base_counts = [{'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'N' : 0, 'DS' : 0} for _ in range(ref_len)]

        # Count the bases (very slowly)
        for pileupcolumn in samfile.pileup(ref, min_base_quality=0, min_mapping_quality=0, max_depth=100000000, stepper="nofilter"):
            ref_pos = pileupcolumn.reference_pos # type: ignore

            for pileupread in pileupcolumn.pileups: # type: ignore
                if not pileupread.is_del and not pileupread.is_refskip:
                    base_counts[ref_pos][pileupread.alignment.query_sequence[pileupread.query_position]] += 1
                else:
                    base_counts[ref_pos]['DS'] += 1

        # Generate per-position statistics
        ref_data = []
        for reference_pos, b_c in enumerate(base_counts):            
            row = []
            coverage = sum(b_c.values())
            if coverage != 0:
                base_probabilities = [count / coverage for count in b_c.values()]
                base_percentages = [100 * probability for probability in base_probabilities]
                entropy = normalising_factor * sum([-(x * math.log2(x)) if x != 0 else 0 for x in base_probabilities])
            else:
                base_percentages = invalid_base_percentages
                entropy = 1

            row.append(reference_pos + 1) # Output one-based coordinates
            row.append(coverage)
            row.extend(b_c.values())
            row.extend(base_percentages)
            row.append(entropy)
            ref_data.append(row)
        data[ref] = ref_data, num_reads

    samfile.close()
    return data


@pytest.mark.parametrize("bam", bams)
def test_basecount(bam):
    # Create an index (if it doesn't already exist) in the same dir as the BAM
    if not os.path.isfile(bam + '.bai'):
        pysam.index(bam) # type: ignore
    
    basecount_data = basecount.get_data(bam)
    test_data = get_test_data(bam)

    # Test for matching references
    assert basecount_data.keys() == test_data.keys()

    # Compare each reference
    for ref in basecount_data.keys():
        basecount_ref_data, basecount_ref_num_reads = basecount_data[ref]
        test_ref_data, test_ref_num_reads = test_data[ref]

        # Test for matching total of reads
        assert basecount_ref_num_reads == test_ref_num_reads

        # Iterate through rows, comparing each
        for basecount_row, test_row in zip(basecount_ref_data, test_ref_data):
            assert basecount_row == test_row
