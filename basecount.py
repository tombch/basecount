import os
import math
import pysam
import argparse
import numpy as np


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


def calculate(bam, references=None):
    data = {}
    normalising_factor = (1 / math.log2(5)) # Normalises maximum entropy to 1
    invalid_base_percentages = [-1, -1, -1, -1, -1] # Percentage values where coverage is zero

    samfile = pysam.AlignmentFile(bam, "rb")

    if references is None:
        references = samfile.references

    for i, ref in enumerate(references):
        ref_data = []

        ref_len = samfile.lengths[i]
        reads = [(read.reference_start, read.cigartuples, read.query_alignment_sequence) for read in samfile.fetch(ref)]

        num_reads = len(reads)

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

        # attempted vectorising
        # base_count_vector = np.array([list(b_c.values()) for b_c in base_counts])
        # coverage_vector = np.sum(base_count_vector, axis=1)
        # base_probabilities_vector = base_count_vector / coverage_vector[:, None]
        # base_percentages_vector = 100 * base_probabilities_vector
        # base_percentages_vector[np.isnan(base_percentages_vector)] = -1

        for reference_pos, b_c in enumerate(base_counts):
            row = []
            coverage = sum(b_c.values()) # type: ignore
            if coverage != 0:
                base_probabilities = [count / coverage for count in b_c.values()]
                base_percentages = [100 * probability for probability in base_probabilities]
                entropy = normalising_factor * sum([-(x * math.log2(x)) if x != 0 else 0 for x in base_probabilities])
                entropy_per_read = entropy / coverage
            else:
                base_percentages = invalid_base_percentages
                entropy = 1
                entropy_per_read = 1

            row.append(reference_pos)
            row.append(coverage)
            row.extend(b_c.values())
            row.extend(base_percentages)
            row.append(entropy)
            row.append(entropy_per_read)
            ref_data.append(row)
        data[ref] = ref_data, num_reads

    samfile.close()
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='path to BAM file')
    parser.add_argument('--references', default=None, action='append', help='reference name(s)')
    parser.add_argument('--bed', default=None, help='path to BED file')
    parser.add_argument('--summarise', default=None, action='store_true', help='display summarising stats instead of per-position stats')
    parser.add_argument('--decimal-places', default=3, type=int, help='default = 3')
    args = parser.parse_args()

    # Create an index (if it doesn't already exist) in the same dir as the BAM
    if not os.path.isfile(args.bam + '.bai'):
        pysam.index(args.bam) # type: ignore

    # Generate data
    data = calculate(args.bam, args.references)

    columns = [
        'reference_position', 
        'num_reads', 
        'num_a', 
        'num_c', 
        'num_g', 
        'num_t', 
        'num_deletions',
        'pc_a', 
        'pc_c', 
        'pc_g', 
        'pc_t', 
        'pc_deletions', 
        'entropy', 
        'entropy_per_read'
    ]
    summary_columns = [
        'ref_length',
        'num_reads',
        'avg_coverage', 
        'num_pos_no_coverage', 
        'pc_pos_no_coverage', 
        'avg_num_deletions', 
        'avg_pc_deletions', 
        'avg_entropy', 
        'avg_entropy_per_read', 
        'mean_entropy_tile_vector',
        'median_entropy_tile_vector',
        'mean_entropy_per_read_tile_vector',
        'median_entropy_per_read_tile_vector'
    ]

    if args.summarise is None:
        # Print all results (in the correct order)
        print('reference_name', '\t'.join(columns), sep='\t')
        for ref, (ref_data, _) in data.items():
            for i, row in enumerate(ref_data):
                print(ref, '\t'.join([str(round(x, args.decimal_places)) for x in row]), sep='\t')
    else:
        # Generate summary statistics
        num_reads_column = columns.index('num_reads')
        num_deletions_column = columns.index('num_deletions')
        percentages_deletions_column = columns.index('pc_deletions')
        entropies_column = columns.index('entropy')
        entropy_per_read_column = columns.index('entropy_per_read')

        for ref, (ref_data, total_reads) in data.items():

            ref_length = len(ref_data)

            num_no_coverage = len([row[num_reads_column] for row in ref_data if row[num_reads_column] == 0])
            pc_no_coverage = 100 * num_no_coverage / ref_length

            num_reads = [row[num_reads_column] for row in ref_data]
            avg_coverage = np.mean(num_reads)

            num_deletions = [row[num_deletions_column] for row in ref_data]
            avg_num_deletions = np.mean(num_deletions)
            
            pc_deletions = [row[percentages_deletions_column] for row in ref_data]
            avg_pc_deletions = np.mean(pc_deletions)

            entropies = [row[entropies_column] for row in ref_data]
            avg_entropies = np.mean(entropies)

            entropy_per_read = [row[entropy_per_read_column] for row in ref_data]
            avg_entropy_per_read = np.mean(entropy_per_read)  
            
            if args.bed is not None:
                scheme = load_scheme(args.bed)
                tile_starts = [tile[2]["inside_start"] for tile in scheme]
                tile_ends = [tile[2]["inside_end"] for tile in scheme]
            else:
                scheme = []
                tile_starts = []
                tile_ends = []

            entropy_data = [[] for _ in scheme]
            entropy_per_read_data = [[] for _ in scheme]
            mean_entropy_vector = []
            median_entropy_vector = []
            mean_entropy_per_read_vector = []
            median_entropy_per_read_vector = []

            for i, (start, end) in enumerate(zip(tile_starts, tile_ends)):
                for j, (entropy, entropy_p_r) in enumerate(zip(entropies, entropy_per_read)):
                    if start <= j <= end:
                        entropy_data[i].append(entropy)
                        entropy_per_read_data[i].append(entropy_p_r)
                
                if entropy_data[i]:
                    mean_entropy_vector.append(np.mean(entropy_data[i]))
                    median_entropy_vector.append(np.median(entropy_data[i]))
                else:
                    mean_entropy_vector.append(-1)
                    median_entropy_vector.append(-1)
                
                if entropy_per_read_data[i]:
                    mean_entropy_per_read_vector.append(np.mean(entropy_per_read_data[i]))
                    median_entropy_per_read_vector.append(np.median(entropy_per_read_data[i]))
                else:
                    mean_entropy_per_read_vector.append(-1)
                    median_entropy_per_read_vector.append(-1)
            
            # Format and display summary statistics
            summary_stats = [ 
                str(round(x, args.decimal_places)) for x in [
                    ref_length, 
                    total_reads,
                    avg_coverage, 
                    num_no_coverage, 
                    pc_no_coverage, 
                    avg_num_deletions, 
                    avg_pc_deletions, 
                    avg_entropies, 
                    avg_entropy_per_read
                ]
            ]

            for vector in [mean_entropy_vector, median_entropy_vector, mean_entropy_per_read_vector, median_entropy_per_read_vector]:
                if vector:
                    summary_stats.append(', '.join([str(round(x, args.decimal_places)) for x in vector]))
                else:
                    summary_stats.append('-')

            print('ref_name', ref, sep='\t')
            for name, val in zip(summary_columns, summary_stats):
                print(name, val, sep='\t')


if __name__ == '__main__':
    main()
