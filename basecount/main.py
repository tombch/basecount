import math
import pysam
import argparse
import numpy as np
from count import bcount


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


# The difference barely matters 
# But still, why compute these ~30000 times when you can just do it once
NORMALISING_FACTOR = 1 / math.log2(5)
SECONDARY_NORMALISING_FACTOR = 1 / math.log2(4)
INVALID_PERCENTAGES = [-1, -1, -1, -1, -1, -1]

def get_position_stats(base_count):
    # Defaults if the coverage is zero
    base_percentages = INVALID_PERCENTAGES
    entropy = 1
    secondary_entropy = 1

    coverage = sum(base_count)
    if coverage != 0:
        base_probabilities = [count / coverage for count in base_count]
        base_percentages = [100 * probability for probability in base_probabilities]
        entropy = NORMALISING_FACTOR * sum([-(x * math.log2(x)) if x != 0 else 0 for x in base_probabilities])

        secondary_base_count = list(base_count)
        secondary_base_count.pop(np.argmax(base_count))
        secondary_coverage = sum(secondary_base_count)
        if secondary_coverage != 0:
            secondary_base_probabilities = [count / secondary_coverage for count in secondary_base_count]
            secondary_entropy = SECONDARY_NORMALISING_FACTOR * sum([-(x * math.log2(x)) if x != 0 else 0 for x in secondary_base_probabilities])
    
    return coverage, base_percentages, entropy, secondary_entropy


def get_data(bam, references=None, long_table=False):
    # Temporarily suppress HTSlib's messages when opening the file
    old_verbosity = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(bam, mode="rb")
    pysam.set_verbosity(old_verbosity) 

    # Iterator through all reads in the file (including unmapped)
    all_reads = samfile.fetch(until_eof=True) 

    if references is None:
        # Calculate data on all references, if no specific ones were given
        references = samfile.references
    else:
        # Determine validity of user-provided references
        for reference in references:
            if not (reference in samfile.references):
                raise Exception(f"{reference} is not a valid reference")
    references = set(references)

    basecount_data = {}
    read_data = {}

    for read in all_reads:
        if not read.is_unmapped:
            ref = read.reference_name
            if ref in references:
                # Create dictionary for a reference if it doesn't already exist
                if read_data.get(ref) is None:
                    read_data[ref] = {}
                    read_data[ref]["ref_len"] = samfile.lengths[samfile.references.index(ref)]
                    read_data[ref]["reads"] = []
                    read_data[ref]["starts"] = []
                    read_data[ref]["ctuples"] = []

                read_data[ref]["reads"].append(read.query_alignment_sequence)
                read_data[ref]["starts"].append(read.reference_start)
                read_data[ref]["ctuples"].append(read.cigartuples)
        
    for ref, ref_data in read_data.items():
        if not (len(ref_data["reads"]) == len(ref_data["starts"]) and len(ref_data["starts"]) == len(ref_data["ctuples"])):
            raise Exception(f"The number of reads, starts and ctuples for {ref} do not match")
        ref_data["num_reads"] = len(ref_data["reads"])

    for ref in references:
        # Count the bases
        base_counts = bcount(
            read_data[ref]["ref_len"], 
            read_data[ref]["reads"], 
            read_data[ref]["starts"], 
            read_data[ref]["ctuples"]
        )

        # Generate per-position statistics
        if not long_table:
            # Wide format (default)
            ref_data = []
            for reference_pos, b_c in enumerate(base_counts):
                coverage, base_percentages, entropy, secondary_entropy = get_position_stats(b_c)
                row = []
                row.append(reference_pos + 1) # Output one-based coordinates
                row.append(coverage)
                row.extend(b_c)
                row.extend(base_percentages)
                row.append(entropy)
                row.append(secondary_entropy)
                ref_data.append(row)
                
            basecount_data[ref] = ref_data, read_data[ref]["num_reads"]
        else:
            # Long format
            ref_data = []
            for reference_pos, b_c in enumerate(base_counts):            
                coverage, base_percentages, entropy, secondary_entropy = get_position_stats(b_c)
                for base, count, percentage in zip(["A", "C", "G", "T", "DS", "N"], b_c, base_percentages):
                    row = []
                    row.append(reference_pos + 1) # Output one-based coordinates
                    row.append(coverage)
                    row.append(base)
                    row.append(count)
                    row.append(percentage)
                    row.append(entropy)
                    row.append(secondary_entropy)
                    ref_data.append(row)

            basecount_data[ref] = ref_data, read_data[ref]["num_reads"]

    samfile.close()
    return basecount_data


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="path to BAM file")
    parser.add_argument("--references", default=None, nargs="+", action="append", help="reference name(s)")
    parser.add_argument("--bed", default=None, nargs="+", action="append", help="path to BED file") 
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--summarise", default=False, action="store_true", help="display summarising stats instead of per-position stats")
    group.add_argument("--long", default=False, action="store_true", help="display per-position stats in long format")
    parser.add_argument("--dp", default=None, nargs="+", action="append", help="default = 3")
    args = parser.parse_args()

    # Move input references into single list, or set as None if none were given
    if args.references is None:
        references = None
    else:
        references = list({ref for ref_list in args.references for ref in ref_list})
    
    # Move bed files into single list, check only one was given
    if args.bed is None:
        bed = None
    else:
        bed_list = list({bed for bed_list in args.bed for bed in bed_list})
        if len(bed_list) > 1:
            raise Exception("Only one BED file can be provided")
        else:
            bed = bed_list[0]
    
    # Handle dp argument
    if args.dp is None:
        dp = 3
    else:
        dp_list = list({dp for dp_list in args.dp for dp in dp_list})
        if len(dp_list) > 1:
            raise Exception("Only one dp value can be provided")
        else:
            dp = int(dp_list[0])

    # Generate data
    data = get_data(args.bam, references=references, long_table=args.long)

    # Output columns
    if not args.long:
        # Wide format columns (default)
        columns = ["position", "coverage", "num_a", "num_c", "num_g", "num_t", "num_ds", "num_n", "pc_a", "pc_c", "pc_g", "pc_t", "pc_ds", "pc_n", "entropy", "secondary_entropy"]
    else:
        # Long format columns (not compatible with --summarise)
        columns = ["position", "coverage", "base", "count", "percentage", "entropy", "secondary_entropy"]

    if not args.summarise:
        if not args.long:
            # Wide format (default)
            print("reference", "\t".join(columns), sep="\t")
            for ref, (ref_data, _) in data.items():
                for i, row in enumerate(ref_data):
                    print(ref, "\t".join([str(round(x, dp)) for x in row]), sep="\t")
        else:
            # Long format
            print("reference", "\t".join(columns), sep="\t")
            for ref, (ref_data, _) in data.items():
                for i, row in enumerate(ref_data):
                    print(ref, "\t".join([str(round(x, dp)) if not isinstance(x, str) else x for x in row]), sep="\t")
    else:
        # Generate summary statistics
        coverage_column = columns.index("coverage")
        entropies_column = columns.index("entropy")

        for ref, (ref_data, total_reads) in data.items():
            ref_length = len(ref_data)
            coverages = [row[coverage_column] for row in ref_data]

            num_ref_coverage = len([cov for cov in coverages if cov != 0])
            pc_ref_coverage = 100 * (num_ref_coverage / ref_length)
            avg_coverage = np.mean([row[coverage_column] for row in ref_data])
            entropies = [row[entropies_column] for row in ref_data]
            avg_entropies = np.mean(entropies)
            
            if bed is not None:
                scheme = load_scheme(bed)
                tile_starts = [tile[2]["inside_start"] for tile in scheme]
                tile_ends = [tile[2]["inside_end"] for tile in scheme]
            else:
                scheme = tile_starts = tile_ends = []

            # Create tile vectors
            coverage_data = [[] for _ in scheme]
            mean_coverage_vector = []
            median_coverage_vector = []
            entropy_data = [[] for _ in scheme]
            mean_entropy_vector = []
            median_entropy_vector = []

            for i, (start, end) in enumerate(zip(tile_starts, tile_ends)):
                for j, (coverage, entropy) in enumerate(zip(coverages, entropies)):
                    if start <= j <= end:
                        coverage_data[i].append(coverage)
                        entropy_data[i].append(entropy)
                
                if coverage_data[i]:
                    mean_coverage_vector.append(np.mean(coverage_data[i]))
                    median_coverage_vector.append(np.median(coverage_data[i]))
                else:
                    mean_coverage_vector.append(-1)
                    median_coverage_vector.append(-1)  
  
                if entropy_data[i]:
                    mean_entropy_vector.append(np.mean(entropy_data[i]))
                    median_entropy_vector.append(np.median(entropy_data[i]))
                else:
                    mean_entropy_vector.append(-1)
                    median_entropy_vector.append(-1)
                
            # Format summary statistics
            summary_stats = {
                "ref_name" : ref,
                "ref_length" : round(ref_length, dp),
                "num_reads" : round(total_reads, dp),
                "pc_ref_coverage" : round(pc_ref_coverage, dp),
                "avg_coverage" : round(avg_coverage, dp),
                "avg_entropy" : round(avg_entropies, dp),
                "mean_coverage_tile_vector" : ", ".join([str(round(x, dp)) for x in mean_coverage_vector]) if mean_coverage_vector else "-",
                "median_coverage_tile_vector" : ", ".join([str(round(x, dp)) for x in median_coverage_vector]) if median_coverage_vector else "-",
                "mean_entropy_tile_vector" : ", ".join([str(round(x, dp)) for x in mean_entropy_vector]) if mean_entropy_vector else "-",
                "median_entropy_tile_vector" : ", ".join([str(round(x, dp)) for x in median_entropy_vector]) if median_entropy_vector else "-"
            }

            # Display summary statistics
            for name, val in summary_stats.items():
                print(name, val, sep="\t")
