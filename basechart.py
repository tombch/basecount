import os
import pysam
import argparse
import concurrent.futures
import math


def calculate_region_percentages(bam, ref, dp, start, end):
    region = []
    samfile = pysam.AlignmentFile(bam, "rb") # type: ignore
    for pileupcolumn in samfile.pileup(ref, start, end, min_base_quality=0):
        if start <= pileupcolumn.pos < end:
            base_count = {
                'A' : 0,
                'C' : 0,
                'G' : 0,
                'T' : 0,
            }
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # TODO: Could there be non-ACGT ?
                    base_count[pileupread.alignment.query_sequence[pileupread.query_position].upper()] += 1

            row = []
            coverage = pileupcolumn.n
            base_probabilities = [count / coverage for count in base_count.values()]
            base_percentages = [round(100 * probability, dp) for probability in base_probabilities]

            row.append(str(pileupcolumn.pos))
            row.extend([str(x) for x in base_percentages])
            row.append(str(round(100 * sum(base_count.values()) / coverage, dp)))
            row.append(str(sum([-(x * math.log2(x)) if x != 0 else 0 for x in base_probabilities])))

            region.append('\t'.join(row))

    samfile.close()
    return '\n'.join(region)


def calculate_region_totals(bam, ref, start, end):
    region = []
    samfile = pysam.AlignmentFile(bam, "rb") # type: ignore
    for pileupcolumn in samfile.pileup(ref, start, end, min_base_quality=0):
        if start <= pileupcolumn.pos < end:
            base_count = {
                'A' : 0,
                'C' : 0,
                'G' : 0,
                'T' : 0,
            }
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # TODO: Could there be non-ACGT ?
                    base_count[pileupread.alignment.query_sequence[pileupread.query_position].upper()] += 1

            row = []

            row.append(str(pileupcolumn.pos))
            row.extend([str(count) for count in base_count.values()])
            row.append(str(sum(base_count.values())))
            
            region.append('\t'.join(row))

    samfile.close()
    return '\n'.join(region)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, help='path to BAM file')
    parser.add_argument('--ref', required=True, help='reference name')
    parser.add_argument('--index', default=None, help='path to index file')
    parser.add_argument('--totals', default=False, action='store_true', help='display base counts instead of percentages')
    parser.add_argument('--decimal-places', default=2, type=int, help='default = 2')
    parser.add_argument('--max-processes', default=None, type=int, help='default = number of processors that the machine has')
    args = parser.parse_args()

    # TODO: Not getting 100% total sum of acgt
    # TODO: Add --summarise

    # Create an index if it doesn't already exist
    if (args.index is None) and (not os.path.isfile(args.bam + '.bai')):
        pysam.index(args.bam) # type: ignore

    # Determine number of bases covered
    # TODO: does pysam have a better way?
    samfile = pysam.AlignmentFile(args.bam, "rb") # type: ignore
    last_pos = [pileupcolumn.pos for pileupcolumn in samfile.pileup(args.ref, min_base_quality=0)][-1]
    samfile.close()

    # Divide genome into appropriate regions to execute in parallel
    region_size = 1000
    regions = [(i, i + region_size) for i in range(0, last_pos, region_size)]

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_processes) as executor:
        if not args.totals:
            # Calculate base percentages at each position
            print('\t'.join(['pos', 'pc_a', 'pc_c', 'pc_g', 'pc_t', 'pc_acgt', 'entropy']))
            results = {start : executor.submit(calculate_region_percentages, args.bam, args.ref, args.decimal_places, start, end) for start, end in regions}
        else:
            # Calculate base totals at each position
            print('\t'.join(['pos', 'num_a', 'num_c', 'num_g', 'num_t', 'num_acgt']))
            results = {start : executor.submit(calculate_region_totals, args.bam, args.ref, start, end) for start, end in regions}

    # Print results in the correct order
    for start, _ in regions:
        print(results[start].result())
    

if __name__ == '__main__':
    main()
