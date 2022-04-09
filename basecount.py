import os
import math
import pysam
import argparse
import statistics


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


def calculate(bam, ref):
    samfile = pysam.AlignmentFile(bam, "rb")
    ref_len = samfile.lengths[0]
    reads = [(read.reference_start, read.cigartuples, read.query_alignment_sequence) for read in samfile.fetch(ref)]
    samfile.close()

    base_counts = [{'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'DEL' : 0} for _ in range(ref_len)]
    operations = {0, 1, 2} # Codes for match, insertion and deletion

    for read in reads:
        start, c_tuples, seq = read
        # Filter out all operations except match, insertion, deletion
        # TODO: keep it simple for now, but this probably has issues
        c_tuples = [c_tuple for c_tuple in c_tuples if c_tuple[0] in operations] # type: ignore

        ref_pos = start
        seq_pos = 0
        for operation, length in c_tuples:
            # Matching
            if operation == 0:
                for r_p, s_p in zip(range(ref_pos, ref_pos + length), range(seq_pos, seq_pos + length)):
                    base_counts[r_p][seq[s_p]] += 1 # type: ignore
                seq_pos += length
                ref_pos += length
            # Insertion in read
            elif operation == 1:
                seq_pos += length
            # Deletion in read
            else:
                for r_p in range(ref_pos, ref_pos + length):
                    base_counts[r_p]['DEL'] += 1 # type: ignore
                ref_pos += length
    
    data = []
    normalising_factor = (1 / math.log2(5)) # Normalises maximum entropy to 1
    for reference_pos, b_c in enumerate(base_counts):
        row = []
        coverage = sum(b_c.values()) # type: ignore
        if coverage != 0:
            base_probabilities = [count / coverage for count in b_c.values()]
            base_percentages = [100 * probability for probability in base_probabilities]
            entropy = normalising_factor * sum([-(x * math.log2(x)) if x != 0 else 0 for x in base_probabilities])
            entropy_per_read = entropy / coverage
        else:
            base_percentages = [-1 for _ in b_c.values()]
            entropy = 1
            entropy_per_read = 1

        row.append(reference_pos)
        row.append(coverage)
        row.extend(b_c.values())
        row.extend(base_percentages)
        row.append(entropy)
        row.append(entropy_per_read)
        data.append(row)

    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, help='path to BAM file')
    parser.add_argument('--ref', required=True, help='reference name')
    parser.add_argument('--bed', default=None, help='path to BED file')
    parser.add_argument('--summarise', default=None, action='store_true', help='display summarising stats instead of per-position stats')
    parser.add_argument('--decimal-places', default=3, type=int, help='default = 3')
    args = parser.parse_args()

    # Create an index (if it doesn't already exist) in the same dir as the BAM
    if not os.path.isfile(args.bam + '.bai'):
        pysam.index(args.bam) # type: ignore

    # Generate data
    data = calculate(args.bam, args.ref)

    columns = ['reference_position', 'num_reads', 'num_a', 'num_c', 'num_g', 'num_t', 'num_deletions','pc_a', 'pc_c', 'pc_g', 'pc_t', 'pc_deletions', 'entropy', 'entropy_per_read']
    summary_columns = ['ref_name', 'ref_length', 'num_no_coverage', 'pc_no_coverage', 'avg_num_reads', 'avg_num_deletions', 'avg_pc_deletions', 'avg_entropy', 'avg_entropy_per_read', 'median_entropy_tile_vector', 'median_entropy_per_read_tile_vector']

    if args.summarise is None:
        # Print all results (in the correct order)
        print('\t'.join(columns))
        for i, row in enumerate(data):
            print('\t'.join([str(round(x, args.decimal_places)) for x in row]))
    else:
        # Generate and display summarising statistics
        if args.bed is not None:
            scheme = load_scheme(args.bed)
            tile_starts = [tile[2]["inside_start"] for tile in scheme]
            tile_ends = [tile[2]["inside_end"] for tile in scheme]
        else:
            scheme = []
            tile_starts = []
            tile_ends = []

        num_reads_column = columns.index('num_reads')
        num_deletions_column = columns.index('num_deletions')
        percentages_deletions_column = columns.index('pc_deletions')
        entropies_column = columns.index('entropy')
        entropy_per_read_column = columns.index('entropy_per_read')

        ref_length = len(data)

        num_no_coverage = len([row[num_reads_column] for row in data if row[num_reads_column] == 0])
        pc_no_coverage = round(100 * num_no_coverage / ref_length, args.decimal_places)

        num_reads = [row[num_reads_column] for row in data]
        avg_num_reads = round(sum(num_reads) / ref_length, args.decimal_places)

        num_deletions = [row[num_deletions_column] for row in data]
        avg_num_deletions = round(sum(num_deletions) / ref_length, args.decimal_places)
        
        pc_deletions = [row[percentages_deletions_column] for row in data]
        avg_pc_deletions = round(sum(pc_deletions) / ref_length, args.decimal_places)

        entropies = [row[entropies_column] for row in data]
        avg_entropies = round(sum(entropies) / ref_length, args.decimal_places)

        entropy_per_read = [row[entropy_per_read_column] for row in data]
        avg_entropy_per_read = round(sum(entropy_per_read) / ref_length, args.decimal_places)  
        
        entropy_data = [[] for _ in scheme]
        entropy_per_read_data = [[] for _ in scheme]
        entropy_vector = []
        entropy_per_read_vector = []
        for i, (start, end) in enumerate(zip(tile_starts, tile_ends)):
            for j, (entropy, entropy_p_r) in enumerate(zip(entropies, entropy_per_read)):
                if start <= j <= end:
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

        summary_stats = [args.ref, ref_length, num_no_coverage, pc_no_coverage, avg_num_reads, avg_num_deletions, avg_pc_deletions, avg_entropies, avg_entropy_per_read]

        for name, val in zip(summary_columns, [str(x) for x in summary_stats] + [entropy_vector_str, entropy_per_read_vector_str]):
            print(name + '\t' + val)


if __name__ == '__main__':
    main()
