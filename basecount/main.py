import math
import pysam
import argparse
import numpy as np
from count import bcount
from basecount.scheme import load_scheme


def get_entropy(probabilities):
    return sum([-(x * math.log2(x)) if x != 0 else 0 for x in probabilities])


def get_basecounts(bam, references=None, long_format=False):
    # Temporarily suppress HTSlib's messages when opening the file
    old_verbosity = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(bam, mode="rb")
    pysam.set_verbosity(old_verbosity) 

    # Determine which references will be covered
    if references is None:
        # Calculate data on all references by default
        references = samfile.references
    else:
        # If references are provided, determine their validity
        for reference in references:
            if not (reference in samfile.references):
                raise Exception(f"{reference} is not a valid reference")
    references = set(references)

    # Iterator through all reads in the file (including unmapped)
    all_reads = samfile.fetch(until_eof=True) 
    
    # Filter for mapped reads only, and store each read's sequence, start and ctuples
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
    
    # For each reference, verify the read data and calculate num_reads
    for ref, ref_data in read_data.items():
        if not (len(ref_data["reads"]) == len(ref_data["starts"]) and len(ref_data["starts"]) == len(ref_data["ctuples"])):
            raise Exception(f"The number of reads, starts and ctuples for {ref} do not match")
        ref_data["num_reads"] = len(ref_data["reads"])

    NORMALISING_FACTOR = 1 / math.log2(6)
    SECONDARY_NORMALISING_FACTOR = 1 / math.log2(5)
    INVALID_PERCENTAGES = [-1, -1, -1, -1, -1, -1]

    basecount_data = {}
    for ref in references:
        # Count the bases
        base_counts = bcount(
            read_data[ref]["ref_len"], 
            read_data[ref]["reads"], 
            read_data[ref]["starts"], 
            read_data[ref]["ctuples"]
        )

        # Generate per-position statistics
        ref_data = []
        for reference_pos, base_count in enumerate(base_counts):
            # Defaults for when the coverage is zero
            base_percentages = INVALID_PERCENTAGES
            entropy = 1
            secondary_entropy = 1
            coverage = sum(base_count)

            if coverage != 0:
                base_probabilities = [count / coverage for count in base_count]
                base_percentages = [100 * probability for probability in base_probabilities]
                entropy = NORMALISING_FACTOR * get_entropy(base_probabilities)

                secondary_base_count = list(base_count)
                secondary_base_count.pop(np.argmax(base_count))
                secondary_coverage = sum(secondary_base_count)
                if secondary_coverage != 0:
                    secondary_base_probabilities = [count / secondary_coverage for count in secondary_base_count]
                    secondary_entropy = SECONDARY_NORMALISING_FACTOR * get_entropy(secondary_base_probabilities)

            # Create row of basecount data, either in long format or wide format 
            # Wide format is default
            if long_format:
                for base, count, percentage in zip(["A", "C", "G", "T", "DS", "N"], base_count, base_percentages):
                    row = []
                    row.append(ref)
                    row.append(reference_pos + 1) # Output one-based coordinates
                    row.append(coverage)
                    row.append(base)
                    row.append(count)
                    row.append(percentage)
                    row.append(entropy)
                    row.append(secondary_entropy)
                    ref_data.append(row)
            else:
                row = []
                row.append(ref)
                row.append(reference_pos + 1) # Output one-based coordinates
                row.append(coverage)
                row.extend(base_count)
                row.extend(base_percentages)
                row.append(entropy)
                row.append(secondary_entropy)
                ref_data.append(row)
            
        basecount_data[ref] = {
            "rows" : ref_data, 
            "num_reads" : read_data[ref]["num_reads"]
        }

    samfile.close()
    return basecount_data


class BaseCount():
    def __init__(self, bam, references=None, long_format=False):
        '''
        Generate and store basecount data.
        
        Arguments:

        * `bam`: path to BAM file to run `basecount` on.
        * `references`: list of references within BAM file to run `basecount` on. If `None`, `basecount` will be run on every reference.
        * `long_format`: `True`/`False` to determine whether basecount data is stored in long or wide format. Is `False` by default.
        '''
        if long_format:
            self.columns = [
                "reference", 
                "position", 
                "coverage", 
                "base", 
                "count", 
                "percentage", 
                "entropy", 
                "secondary_entropy"
            ]
        else:
            self.columns = [
                "reference", 
                "position", 
                "coverage", 
                "num_a", 
                "num_c", 
                "num_g", 
                "num_t", 
                "num_ds", 
                "num_n", 
                "pc_a", 
                "pc_c", 
                "pc_g", 
                "pc_t", 
                "pc_ds", 
                "pc_n", 
                "entropy", 
                "secondary_entropy"
            ]
        self.data = get_basecounts(bam, references, long_format)
        self.references = list(self.data.keys())
        self.reference_lengths = {ref : len(self.data[ref]["rows"]) for ref in self.references}

    def rows(self, reference=None):
        '''
        Returns an iterator of lists for the basecount data.

        If a reference is given, then only basecount data for that reference will be yielded.
        '''
        if reference is None:
            for ref in self.data.keys():
                for row in self.data[ref]["rows"]:
                    yield row
        else:
            if self.data.get(reference) is None:
                raise Exception(f"{reference} is not a valid reference")
            for row in self.data[reference]["rows"]:
                yield row

    def records(self, reference=None):
        '''
        Returns an iterator of dictionaries (with column names as the keys) for the basecount data.

        If a reference is given, then only basecount data for that reference will be yielded.
        '''
        if reference is None:
            for ref in self.data.keys():
                for row in self.data[ref]["rows"]:
                    yield dict(zip(self.columns, row))
        else:
            if self.data.get(reference) is None:
                raise Exception(f"{reference} is not a valid reference")
            for row in self.data[reference]["rows"]:
                yield dict(zip(self.columns, row))

    def num_reads(self, reference=None):
        '''
        Returns the total number of reads, across all references.

        If a reference is given, returns the total number of reads just for that reference.
        '''
        if reference is None:
            return sum([self.data[ref]["num_reads"] for ref in self.references])
        else:
            if self.data.get(reference) is None:
                raise Exception(f"{reference} is not a valid reference")
            return self.data[reference]["num_reads"]


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

    # Generate the data
    bc = BaseCount(args.bam, references=references, long_format=args.long)

    if not args.summarise:
        print("\t".join(bc.columns))
        for row in bc.rows():
            print("\t".join([str(round(x, dp)) if not isinstance(x, str) else x for x in row]), sep="\t")
    else:
        # Generate summary statistics        
        for ref in bc.references:
            coverages = []
            entropies = []
            secondary_entropies = []

            for record in bc.records(ref):
                coverages.append(record["coverage"])
                entropies.append(record["entropy"])
                secondary_entropies.append(record["secondary_entropy"])

            avg_coverage = np.mean(coverages)
            avg_entropies = np.mean(entropies)

            ref_length = bc.reference_lengths[ref]
            pc_ref_coverage = 100 * (len([cov for cov in coverages if cov != 0]) / ref_length)
            
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
            
            secondary_entropy_data = [[] for _ in scheme]
            mean_secondary_entropy_vector = []
            median_secondary_entropy_vector = []

            for i, (start, end) in enumerate(zip(tile_starts, tile_ends)):
                for j, (coverage, entropy, secondary_entropy) in enumerate(zip(coverages, entropies, secondary_entropies)):
                    if start <= j <= end:
                        coverage_data[i].append(coverage)
                        entropy_data[i].append(entropy)
                        secondary_entropy_data[i].append(secondary_entropy)
                
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

                if secondary_entropy_data[i]:
                    mean_secondary_entropy_vector.append(np.mean(secondary_entropy_data[i]))
                    median_secondary_entropy_vector.append(np.median(secondary_entropy_data[i]))
                else:
                    mean_secondary_entropy_vector.append(-1)
                    median_secondary_entropy_vector.append(-1)
                
            # Format summary statistics
            summary_stats = {
                "reference" : ref,
                "reference_length" : round(ref_length, dp),
                "num_reads" : round(bc.num_reads(ref), dp),
                "pc_reference_coverage" : round(pc_ref_coverage, dp),
                "avg_depth" : round(avg_coverage, dp),
                "avg_entropy" : round(avg_entropies, dp),
                "mean_coverage_tile_vector" : ", ".join([str(round(x, dp)) for x in mean_coverage_vector]) if mean_coverage_vector else "-",
                "median_coverage_tile_vector" : ", ".join([str(round(x, dp)) for x in median_coverage_vector]) if median_coverage_vector else "-",
                "mean_entropy_tile_vector" : ", ".join([str(round(x, dp)) for x in mean_entropy_vector]) if mean_entropy_vector else "-",
                "median_entropy_tile_vector" : ", ".join([str(round(x, dp)) for x in median_entropy_vector]) if median_entropy_vector else "-",
                "mean_secondary_entropy_tile_vector" : ", ".join([str(round(x, dp)) for x in mean_secondary_entropy_vector]) if mean_secondary_entropy_vector else "-",
                "median_secondary_entropy_tile_vector" : ", ".join([str(round(x, dp)) for x in median_secondary_entropy_vector]) if median_secondary_entropy_vector else "-"
            }

            # Display summary statistics
            for name, val in summary_stats.items():
                print(name, val, sep="\t")
