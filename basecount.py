import os
import math
import pysam
import argparse
import numpy as np
import basecomp


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
                    # Push the window region to the leftmost left position
                    tiles_dict[tile]["start"] = start
                if end > tiles_dict[tile]["inside_start"]:
                    # Open the start of the inner window to the rightmost left position
                    tiles_dict[tile]["inside_start"] = end

            elif "RIGHT" in side.upper():
                if tiles_dict[tile]["end"] == -1:
                    tiles_dict[tile]["end"] = end
                    tiles_dict[tile]["inside_end"] = start

                if end > tiles_dict[tile]["end"]:
                    # Stretch the window out to the rightmost right position
                    tiles_dict[tile]["end"] = end
                if start < tiles_dict[tile]["inside_end"]:
                    # Close the end of the inner window to the leftmost right position
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

        tiles_list = sorted(tiles_list, key=lambda x: int(x[1])) # Sort by tile number
        if clip: # Default
            # Iterate through tiles and clip
            new_tiles = []
            for tile_index, tile_tuple in enumerate(tiles_list):
                tile_dict = dict(tile_tuple[2])

                # Clip the start of this window to the end of the last window, if there is a last window
                if tile_index > 0:
                    tile_dict["inside_start"] = tiles_list[tile_index - 1][2]["end"]

                # Clip the end of this window to the start of the next window, if there is a next window
                if tile_index < len(tiles_list) - 1:
                    tile_dict["inside_end"] = tiles_list[tile_index + 1][2]["start"]

                new_tiles.append((tile_tuple[0], tile_tuple[1], tile_dict))
        else:
            new_tiles = tiles_list

    return new_tiles


def get_data(bam, references=None):    
    samfile = pysam.AlignmentFile(bam, "rb")

    # If no references are specified, calculate stats for every contig
    if references is None:
        references = samfile.references

    data = {}
    normalising_factor = 1 / math.log2(5) # Normalises maximum entropy to 1
    invalid_base_percentages = [-1, -1, -1, -1, -1, -1] # Percentage values where coverage is zero

    for ref in references:
        # Prepare read data
        ref_len = samfile.lengths[samfile.references.index(ref)]
        read_data = samfile.fetch(ref) # Iterator through all mapped reads in the contig
        reads = []
        starts = []
        ctuples = []
        for r_d in read_data:
            reads.append(r_d.query_alignment_sequence)
            starts.append(r_d.reference_start)
            ctuples.append(r_d.cigartuples)
        num_reads = len(reads)

        # Count the bases
        base_counts = basecomp.compute(ref_len, reads, starts, ctuples)

        # Generate per-position statistics
        ref_data = []
        for reference_pos, b_c in enumerate(base_counts):
            # TODO: Is it best to ignore these completely?
            # b_c.pop('N')
            
            row = []
            coverage = sum(b_c)
            if coverage != 0:
                base_probabilities = [count / coverage for count in b_c]
                base_percentages = [100 * probability for probability in base_probabilities]
                entropy = normalising_factor * sum([-(x * math.log2(x)) if x != 0 else 0 for x in base_probabilities])
            else:
                base_percentages = invalid_base_percentages
                entropy = 1 # TODO: is this an appropriate value?

            row.append(reference_pos + 1) # Output one-based coordinates
            row.append(coverage)
            row.extend(b_c)
            row.extend(base_percentages)
            row.append(entropy)
            ref_data.append(row)
        data[ref] = ref_data, num_reads

    samfile.close()
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='path to BAM file')
    parser.add_argument('--references', default=None, nargs='+', action='append', help='reference name(s)')
    parser.add_argument('--bed', default=None, nargs='+', action='append', help='path to BED file')
    parser.add_argument('--summarise', default=False, action='store_true', help='display summarising stats instead of per-position stats')
    parser.add_argument('--dp', default=None, nargs='+', action='append', help='default = 3')
    args = parser.parse_args()

    # Move input references into single list, or set as None if none were given
    if args.references is not None:
        references = list({ref for ref_list in args.references for ref in ref_list})
    else:
        references = None
    
    # Move bed files into single list, check only one was given
    if args.bed is not None:
        bed_list = list({bed for bed_list in args.bed for bed in bed_list})
        if len(bed_list) > 1:
            raise Exception('Only one BED file can be provided')
        else:
            bed = bed_list[0]
    else:
        bed = None
    
    # Handle dp argument
    if args.dp is None:
        dp = 3
    else:
        dp_list = list({dp for dp_list in args.dp for dp in dp_list})
        if len(dp_list) > 1:
            raise Exception('Only one dp value can be provided')
        else:
            dp = int(dp_list[0])

    # Create an index (if it doesn't already exist) in the same dir as the BAM
    if not os.path.isfile(args.bam + '.bai'):
        pysam.index(args.bam) # type: ignore

    # Generate data
    data = get_data(args.bam, references)

    # Output columns
    columns = [
        'reference_position', 
        'num_reads', 
        'num_a', 
        'num_c', 
        'num_g', 
        'num_t',
        'num_n',
        'num_deletions_skips',
        'pc_a', 
        'pc_c', 
        'pc_g',
        'pc_t',
        'pc_n',
        'pc_deletions_skips', 
        'entropy', 
    ]

    if not args.summarise:
        # Print all results (in the correct order)
        print('reference_name', '\t'.join(columns), sep='\t')
        for ref, (ref_data, _) in data.items():
            for i, row in enumerate(ref_data):
                print(ref, '\t'.join([str(round(x, dp)) for x in row]), sep='\t')
    else:
        # Generate summary statistics
        num_reads_column = columns.index('num_reads')
        num_deletions_column = columns.index('num_deletions_skips')
        percentages_deletions_column = columns.index('pc_deletions_skips')
        entropies_column = columns.index('entropy')

        for ref, (ref_data, total_reads) in data.items():
            ref_length = len(ref_data)

            num_reads = [row[num_reads_column] for row in ref_data]
            avg_coverage = np.mean(num_reads)

            num_no_coverage = len([row[num_reads_column] for row in ref_data if row[num_reads_column] == 0])
            pc_coverage = 100 - (100 * num_no_coverage / ref_length)

            num_deletions_skips = [row[num_deletions_column] for row in ref_data]
            avg_num_deletions_skips = np.mean(num_deletions_skips)
            
            pc_deletions_skips = [row[percentages_deletions_column] for row in ref_data]
            avg_pc_deletions_skips = np.mean(pc_deletions_skips)

            entropies = [row[entropies_column] for row in ref_data]
            avg_entropies = np.mean(entropies)
            
            if bed is not None:
                scheme = load_scheme(bed)
                tile_starts = [tile[2]["inside_start"] for tile in scheme]
                tile_ends = [tile[2]["inside_end"] for tile in scheme]
            else:
                scheme = []
                tile_starts = []
                tile_ends = []

            # Create tile vectors
            entropy_data = [[] for _ in scheme]
            mean_entropy_vector = []
            median_entropy_vector = []

            for i, (start, end) in enumerate(zip(tile_starts, tile_ends)):
                for j, entropy in enumerate(entropies):
                    if start <= j <= end:
                        entropy_data[i].append(entropy)
                
                if entropy_data[i]:
                    mean_entropy_vector.append(np.mean(entropy_data[i]))
                    median_entropy_vector.append(np.median(entropy_data[i]))
                else:
                    mean_entropy_vector.append(-1)
                    median_entropy_vector.append(-1)
                
            # Format summary statistics
            summary_stats = {
                'ref_name' : ref,
                'ref_length' : round(ref_length, dp),
                'num_reads' : round(total_reads, dp),
                'avg_coverage' : round(avg_coverage, dp),
                'pc_ref_coverage' : round(pc_coverage, dp),
                'num_pos_no_coverage' : round(num_no_coverage, dp), 
                'avg_num_deletions_skips' : round(avg_num_deletions_skips, dp), 
                'avg_pc_deletions_skips' : round(avg_pc_deletions_skips, dp), 
                'avg_entropy' : round(avg_entropies, dp), 
                'mean_entropy_tile_vector' : ', '.join([str(round(x, dp)) for x in mean_entropy_vector]) if mean_entropy_vector else '-',
                'median_entropy_tile_vector' : ', '.join([str(round(x, dp)) for x in median_entropy_vector]) if median_entropy_vector else '-',
            }

            # Display summary statistics
            for name, val in summary_stats.items():
                print(name, val, sep='\t')


if __name__ == '__main__':
    main()
