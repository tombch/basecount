import os
import math
import pysam
import argparse
import statistics
import concurrent.futures


# Written by Sam Nicholls as part of swell
# https://github.com/SamStudio8/swell
def load_scheme(bed, clip=True):
    tiles_dict = {}
    with open(bed) as scheme_fh:
        for line in scheme_fh:
            data = line.strip().split()
            start, end, tile = int(data[1]), int(data[2]), data[3] 
            scheme, tile, side = tile.split("_", 2)

            if tile not in tiles_dict:
                tiles_dict[tile] = {
                    "start": -1,
                    "inside_start": -1,
                    "inside_end": -1,
                    "end": -1,
                }

            if "LEFT" in side.upper():
                if tiles_dict[tile]["start"] == -1:
                    tiles_dict[tile]["start"] = start
                    tiles_dict[tile]["inside_start"] = end

                if start < tiles_dict[tile]["start"]:
                    # push the window region to the leftmost left position
                    tiles_dict[tile]["start"] = start
                if end > tiles_dict[tile]["inside_start"]:
                    # open the start of the inner window to the rightmost left position
                    tiles_dict[tile]["inside_start"] = end

            elif "RIGHT" in side.upper():
                if tiles_dict[tile]["end"] == -1:
                    tiles_dict[tile]["end"] = end
                    tiles_dict[tile]["inside_end"] = start

                if end > tiles_dict[tile]["end"]:
                    # stretch the window out to the rightmost right position
                    tiles_dict[tile]["end"] = end
                if start < tiles_dict[tile]["inside_end"]:
                    # close the end of the inner window to the leftmost right position
                    tiles_dict[tile]["inside_end"] = start
        
        tiles_list = []
        tiles_seen = set()
        scheme_fh.seek(0)
        for line in scheme_fh:
            data = line.strip().split()
            start, end, tile = data[1], data[2], data[3] 
            scheme, tile, side = tile.split("_", 2)
            tile_tup = (scheme, tile, tiles_dict[tile])
            if tiles_dict[tile]["inside_start"] != -1 and tiles_dict[tile]["inside_end"] != -1 and tile not in tiles_seen:
                tiles_list.append(tile_tup)
                tiles_seen.add(tile)

        tiles_list = sorted(tiles_list, key=lambda x: int(x[1])) # sort by tile number
        if clip: # Default
            # iterate through tiles and clip
            new_tiles = []
            for tile_index, tile_tuple in enumerate(tiles_list):
                tile_dict = dict(tile_tuple[2]) # copy the dict god this is gross stuff

                # Clip the start of this window to the end of the last window
                # (if there is a last window)
                if tile_index > 0:
                    tile_dict["inside_start"] = tiles_list[tile_index - 1][2]["end"]

                # Clip the end of this window to the start of the next window
                # (if there is a next window)
                if tile_index < len(tiles_list) - 1:
                    tile_dict["inside_end"] = tiles_list[tile_index + 1][2]["start"]

                new_tiles.append((tile_tuple[0], tile_tuple[1], tile_dict))
        else:
            new_tiles = tiles_list

    return new_tiles


def calculate_region(bam, ref, start, end):
    region = []
    samfile = pysam.AlignmentFile(bam, "rb") # type: ignore
    normalising_factor = (1 / math.log2(5)) # Normalises maximum entropy to 1
    for pileupcolumn in samfile.pileup(ref, start, end, min_base_quality=0):
        reference_pos = pileupcolumn.reference_pos # type: ignore
        if start <= reference_pos < end:
            base_count = {
                'A' : 0,
                'C' : 0,
                'G' : 0,
                'T' : 0,
                'SD' : 0 # skips and deletions
            }
            for pileupread in pileupcolumn.pileups: # type: ignore
                if not pileupread.is_del and not pileupread.is_refskip:
                    # TODO: Could there be non-ACGT ?
                    base_count[pileupread.alignment.query_sequence[pileupread.query_position].upper()] += 1
                else:
                    base_count['SD'] += 1

            row = []
            coverage = pileupcolumn.nsegments # type: ignore
            base_probabilities = [count / coverage for count in base_count.values()]
            base_percentages = [100 * probability for probability in base_probabilities]
            entropy = normalising_factor * sum([-(x * math.log2(x)) if x != 0 else 0 for x in base_probabilities])
            entropy_per_read = entropy / coverage

            row.append(reference_pos)
            row.append(coverage)
            row.extend(base_count.values())
            row.extend(base_percentages)
            row.append(entropy)
            row.append(entropy_per_read)
            region.append(row)

    samfile.close()
    return region


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, help='path to BAM file')
    parser.add_argument('--ref', required=True, help='reference name')
    parser.add_argument('--bed', default=None, help='path to BED file') # TODO
    parser.add_argument('--index', default=None, help='path to index file')
    parser.add_argument('--summarise', default=None, action='store_true', help='display summarising stats instead of per-position stats')
    parser.add_argument('--decimal-places', default=3, type=int, help='default = 3')
    parser.add_argument('--max-processes', default=None, type=int, help='default = number of processors that the machine has')
    args = parser.parse_args()

    # Create an index if it doesn't already exist
    if (args.index is None) and (not os.path.isfile(args.bam + '.bai')):
        pysam.index(args.bam) # type: ignore

    # Determine number of bases covered
    # TODO: does pysam have a better way?
    samfile = pysam.AlignmentFile(args.bam, "rb") # type: ignore
    positions = [pileupcolumn.reference_pos for pileupcolumn in samfile.pileup(args.ref, min_base_quality=0)] # type: ignore
    first_pos = positions[0]
    last_pos = positions[-1]
    samfile.close()

    # Divide genome into appropriate regions to execute in parallel
    if args.max_processes is None:
        region_size = int(last_pos / os.cpu_count())
    else:
        region_size = int(last_pos / args.max_processes)
    regions = [(i, i + region_size) for i in range(0, last_pos + 1, region_size)]
    last_region_end = regions[-1][1]
    if last_region_end > last_pos + 1:
        regions.pop()
        regions[-1] = regions[-1][0], last_pos + 1
    
    # Generate results
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_processes) as executor:
        results = {start : executor.submit(calculate_region, args.bam, args.ref, start, end) for start, end in regions}

    columns = ['reference_pos', 'num_reads', 'num_a', 'num_c', 'num_g', 'num_t', 'num_skip_del','pc_a', 'pc_c', 'pc_g', 'pc_t', 'pc_skip_del', 'entropy', 'entropy_per_read']
    summary_columns = ['ref_name', 'first_ref_pos', 'last_ref_pos', 'avg_num_reads', 'avg_num_skip_del', 'avg_pc_skip_del', 'avg_entropy', 'avg_entropy_per_read', 'median_entropy_tile_vector', 'median_entropy_per_read_tile_vector']

    if args.summarise is None:
        # Print all results (in the correct order)
        print('\t'.join(columns))
        for start, _ in regions:
            print('\n'.join(['\t'.join([str(round(x, args.decimal_places)) for x in result]) for result in results[start].result()]))
    else:
        if args.bed is not None:
            scheme = load_scheme(args.bed)
            tile_starts = [tile[2]["inside_start"] for tile in scheme]
            tile_ends = [tile[2]["inside_end"] for tile in scheme]
        else:
            scheme = []
            tile_starts = []
            tile_ends = []

        # Summarise results
        # TODO: Averages currently taken over the number of positions covered by reads, NOT the length of the reference
        num_reads_column = columns.index('num_reads')
        skips_deletions_column = columns.index('num_skip_del')
        percentages_skips_deletions_column = columns.index('pc_skip_del')
        entropies_column = columns.index('entropy')
        entropy_per_read_column = columns.index('entropy_per_read')
        
        num_reads = [result[num_reads_column] for start, _ in regions for result in results[start].result()]
        avg_num_reads = round(sum(num_reads) / len(num_reads), args.decimal_places)

        skips_deletions = [result[skips_deletions_column] for start, _ in regions for result in results[start].result()]
        avg_skips_deletions = round(sum(skips_deletions) / len(skips_deletions), args.decimal_places)
        
        pc_skips_deletions = [result[percentages_skips_deletions_column] for start, _ in regions for result in results[start].result()]
        avg_pc_skips_deletions = round(sum(pc_skips_deletions) / len(pc_skips_deletions), args.decimal_places)

        entropies = [result[entropies_column] for start, _ in regions for result in results[start].result()]
        avg_entropies = round(sum(entropies) / len(entropies), args.decimal_places)

        entropy_per_read = [result[entropy_per_read_column] for start, _ in regions for result in results[start].result()]
        avg_entropy_per_read = round(sum(entropy_per_read) / len(entropy_per_read), args.decimal_places)  
        
        entropy_data = [[] for _ in scheme]
        entropy_per_read_data = [[] for _ in scheme]
        entropy_vector = []
        entropy_per_read_vector = []
        for i, (start, end) in enumerate(zip(tile_starts, tile_ends)):
            for j, (entropy, entropy_p_r) in enumerate(zip(entropies, entropy_per_read)):
                if start <= first_pos + j <= end:
                    entropy_data[i].append(entropy)
                    entropy_per_read_data[i].append(entropy_p_r)
            
            if entropy_data[i]:
                entropy_vector.append(round(statistics.median(entropy_data[i]), args.decimal_places))
            else:
                entropy_vector.append(-1)
            
            if entropy_per_read_data[i]:
                entropy_per_read_vector.append(round(statistics.median(entropy_per_read_data[i]), args.decimal_places))
            else:
                entropy_per_read_vector.append(-1)
        
        if entropy_vector:
            entropy_vector_str = ','.join([str(x) for x in entropy_vector])
        else:
            entropy_vector_str = '-'

        if entropy_per_read_vector:
            entropy_per_read_vector_str = ','.join([str(x) for x in entropy_per_read_vector])
        else:
            entropy_per_read_vector_str = '-'

        summary_stats = [args.ref, first_pos, last_pos, avg_num_reads, avg_skips_deletions, avg_pc_skips_deletions, avg_entropies, avg_entropy_per_read]

        for name, val in zip(summary_columns, [str(x) for x in summary_stats] + [entropy_vector_str, entropy_per_read_vector_str]):
            print(name + '\t' + val)


if __name__ == '__main__':
    main()
