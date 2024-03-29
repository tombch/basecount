import math
import pysam
import argparse
import numpy as np
from count import bcount
from basecount.scheme import load_scheme
from basecount.version import __version__


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
                secondary_base_probabilities = [
                    count / secondary_coverage for count in secondary_base_count
                ]
                secondary_entropy = secondary_normalising_factor * get_entropy(
                    secondary_base_probabilities
                )

        # Create row of basecount data, either in long format or wide format
        # Wide format is default
        if long_format:
            for base, count, percentage in zip(bases, base_count, base_percentages):
                row = []
                row.append(ref)
                row.append(reference_pos + 1)  # Output one-based coordinates
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
            row.append(reference_pos + 1)  # Output one-based coordinates
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


def init_read_chunk(references):
    read_chunk = {}
    for ref in references:
        read_chunk[ref] = {"reads": [], "qualities": [], "starts": [], "ctuples": []}
    return read_chunk


def get_basecounts(
    bam,
    references=None,
    min_base_quality=0,
    min_mapping_quality=0,
    chunk_size=1000000,
    show_n_bases=False,
    long_format=False,
):
    samfile = open_samfile(bam)
    references = get_references(samfile, references)
    reference_lengths = {
        ref: samfile.lengths[samfile.references.index(ref)] for ref in references
    }

    # Iterator through all reads in the file (including unmapped)
    # This is the only way to access the reads without requiring an index file
    reads = samfile.fetch(until_eof=True)

    basecount_data = {}
    for ref in references:
        basecount_data[ref] = {
            "data": np.zeros((reference_lengths[ref], 6), dtype=int),
            "num_reads": 0,
        }

    # Create empty read chunk
    read_chunk = init_read_chunk(references)
    num_chunk_reads = 0

    # Iterate through all reads
    for read in reads:
        if num_chunk_reads == chunk_size:
            # Process the chunk
            for ref in references:
                # Basecount the reads for the reference
                bcounts = bcount(
                    reference_lengths[ref],
                    min_base_quality,
                    read_chunk[ref]["reads"],
                    read_chunk[ref]["qualities"],
                    read_chunk[ref]["starts"],
                    read_chunk[ref]["ctuples"],
                )
                # Add the data to the current basecount data for the reference
                basecount_data[ref]["data"] = np.add(
                    basecount_data[ref]["data"], bcounts
                )
                basecount_data[ref]["num_reads"] += len(read_chunk[ref]["reads"])

            # Empty the read chunk
            read_chunk = init_read_chunk(references)
            num_chunk_reads = 0

        # Append read to the chunk, if it is mapped and gte to min mapping quality
        if (not read.is_unmapped) and (read.mapping_quality >= min_mapping_quality):
            read_chunk[read.reference_name]["reads"].append(
                read.query_alignment_sequence
            )
            read_chunk[read.reference_name]["qualities"].append(
                read.query_alignment_qualities
            )
            read_chunk[read.reference_name]["starts"].append(read.reference_start)
            read_chunk[read.reference_name]["ctuples"].append(read.cigartuples)
            num_chunk_reads += 1

    # Count the remaining bases left in the chunk
    for ref in references:
        # Basecount the reads for the reference
        bcounts = bcount(
            reference_lengths[ref],
            min_base_quality,
            read_chunk[ref]["reads"],
            read_chunk[ref]["qualities"],
            read_chunk[ref]["starts"],
            read_chunk[ref]["ctuples"],
        )
        # Add the data to the current basecount data for the reference
        basecount_data[ref]["data"] = np.add(basecount_data[ref]["data"], bcounts)
        basecount_data[ref]["num_reads"] += len(read_chunk[ref]["reads"])

    # Generate statistics
    basecount_stats = {}
    for ref in references:
        basecount_stats[ref] = {
            "rows": get_stats(
                basecount_data[ref]["data"].tolist(),
                ref,
                show_n_bases=show_n_bases,
                long_format=long_format,
            ),
            "num_reads": basecount_data[ref]["num_reads"],
        }

    samfile.close()
    return basecount_stats


class BaseCount:
    def __init__(
        self,
        bam,
        references=None,
        min_base_quality=0,
        min_mapping_quality=0,
        chunk_size=1000000,
        show_n_bases=False,
        long_format=False,
    ):
        """
        Generate and store basecount data.

        Arguments:

        * `bam`: Path to BAM file (an index file is not required).
        * `references`: List of specific references to run basecount on. If `None`, `basecount` will be run on every reference.
        * `min_base_quality`: Minimum quality of a base, for the base to be included in the `basecount` data. Default: `0`.
        * `min_mapping_quality`: Minimum quality of a read mapping, for the read to be included in the `basecount` data. Default: `0`.
        * `chunk_size`: Max number of reads loaded into memory and basecounted at a given time. Default: `1000000`.
        * `show_n_bases`: Show counts of `N` bases from reads, and include them in statistics. Default: `False`.
        * `long_format`: `True`/`False` to determine whether `basecount` data is stored in long or wide format. Default: `False`.
        """
        if long_format:
            self.columns = [
                "reference",
                "position",
                "coverage",
                "base",
                "count",
                "percentage",
                "entropy",
                "secondary_entropy",
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
                "secondary_entropy",
            ]
            if not show_n_bases:
                self.columns.pop(self.columns.index("num_n"))
                self.columns.pop(self.columns.index("pc_n"))
        self.data = get_basecounts(
            bam,
            references=references,
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
            chunk_size=chunk_size,
            show_n_bases=show_n_bases,
            long_format=long_format,
        )
        self.references = list(self.data.keys())
        self.reference_lengths = {
            ref: len(self.data[ref]["rows"]) for ref in self.references
        }

    def rows(self, reference=None):
        """
        Returns an iterator of lists for the basecount data.

        If a reference is given, then only basecount data for that reference will be yielded.
        """
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
        """
        Returns an iterator of dictionaries (with column names as the keys) for the basecount data.

        If a reference is given, then only basecount data for that reference will be yielded.
        """
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
        """
        Returns the total number of reads, across all references.

        If a `reference` is given, returns the total number of reads just for that reference.
        """
        if reference is None:
            return sum([self.data[ref]["num_reads"] for ref in self.references])
        else:
            if self.data.get(reference) is None:
                raise Exception(f"{reference} is not a valid reference")
            return self.data[reference]["num_reads"]

    def mean_coverage(self, reference=None):
        """
        Returns the mean coverage across all positions, across all references.

        If a `reference` is given, returns the average coverage across all positions of that reference only.
        """
        coverages = []
        if reference is None:
            for record in self.records():
                coverages.append(record["coverage"])
        else:
            if self.data.get(reference) is None:
                raise Exception(f"{reference} is not a valid reference")
            for record in self.records(reference=reference):
                coverages.append(record["coverage"])
        return np.mean(coverages)

    def mean_entropy(self, reference=None, min_coverage=0):
        """
        Returns the mean entropy across all positions, across all references.

        If a `reference` is given, returns the average entropy across all positions of that reference only.
        """
        entropies = []
        if reference is None:
            for record in self.records():
                if record["coverage"] >= min_coverage:
                    entropies.append(record["entropy"])
        else:
            if self.data.get(reference) is None:
                raise Exception(f"{reference} is not a valid reference")
            for record in self.records(reference=reference):
                if record["coverage"] >= min_coverage:
                    entropies.append(record["entropy"])
        return np.mean(entropies)


def handle_arg(arg, name, default=None, provided_once=False):
    if arg is None:
        # If nothing was provided, return the default value
        value = default
    else:
        if provided_once:
            if len(arg) > 1:
                raise Exception(f"Argument --{name} can only be provided once")
            else:
                value = arg[0]
        else:
            # Gather all values into single list
            value = list({a for a_list in arg for a in a_list})
    return value


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="Path to BAM file (an index file is not required)")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument(
        "--references",
        default=None,
        nargs="+",
        action="append",
        help="Choose specific reference(s) to run basecount on",
    )
    parser.add_argument(
        "--min-base-quality", default=None, action="append", help="Default value: 0"
    )
    parser.add_argument(
        "--min-mapping-quality", default=None, action="append", help="Default value: 0"
    )
    parser.add_argument(
        "--chunk-size",
        default=None,
        action="append",
        help="Max number of reads loaded into memory and basecounted at a given time. Default value: 1000000",
    )
    parser.add_argument(
        "--show-n-bases",
        default=False,
        action="store_true",
        help="Show counts of 'N' bases from reads, and include them in statistics",
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--long-format",
        default=False,
        action="store_true",
        help="Output per-position statistics in long format, instead of the default wide format",
    )
    group.add_argument(
        "--summarise",
        default=False,
        action="store_true",
        help="Output summary statistics",
    )
    group.add_argument(
        "--summarise-with-bed",
        default=None,
        action="append",
        metavar="BED_FILE",
        help="Output summary statistics and amplicon vectors (calculated using the provided BED file)",
    )

    parser.add_argument(
        "--decimal-places", default=None, action="append", help="Default value: 3"
    )

    args = parser.parse_args()

    # Handle arguments
    references = handle_arg(args.references, "references")
    min_base_quality = int(handle_arg(args.min_base_quality, "min-base-quality", default=0, provided_once=True))  # type: ignore
    min_mapping_quality = int(handle_arg(args.min_mapping_quality, "min-mapping-quality", default=0, provided_once=True))  # type: ignore
    chunk_size = int(handle_arg(args.chunk_size, "chunk-size", default=1000000, provided_once=True))  # type: ignore
    bed = handle_arg(args.summarise_with_bed, "bed", provided_once=True)
    decimal_places = int(handle_arg(args.decimal_places, "decimal_places", default=3, provided_once=True))  # type: ignore

    # Generate the data
    bc = BaseCount(
        args.bam,
        references=references,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        chunk_size=chunk_size,
        show_n_bases=args.show_n_bases,
        long_format=args.long_format,
    )

    if (not args.summarise) and (bed is None):
        # Display per-position statistics
        print("\t".join(bc.columns))
        for row in bc.rows():
            print(
                "\t".join(
                    [
                        str(round(x, decimal_places)) if not isinstance(x, str) else x
                        for x in row
                    ]
                ),
                sep="\t",
            )
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
            pc_ref_coverage = 100 * (
                len([cov for cov in coverages if cov != 0]) / ref_length
            )

            # Format summary statistics
            summary_stats = {
                "reference_name": ref,
                "reference_length": round(ref_length, decimal_places),
                "num_reads": round(bc.num_reads(ref), decimal_places),
                "pc_reference_coverage": round(pc_ref_coverage, decimal_places),
                "avg_depth": round(avg_coverage, decimal_places),
                "avg_entropy": round(avg_entropies, decimal_places),
            }

            # Display summary statistics
            for name, val in summary_stats.items():
                print(name, val, sep="\t")

            if bed is not None:
                scheme = load_scheme(bed)
                tile_starts = [tile[2]["inside_start"] for tile in scheme]
                tile_ends = [tile[2]["inside_end"] for tile in scheme]

                # Generate amplicon vectors
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
                    for j, (coverage, entropy, secondary_entropy) in enumerate(
                        zip(coverages, entropies, secondary_entropies)
                    ):
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
                        mean_secondary_entropy_vector.append(
                            np.mean(secondary_entropy_data[i])
                        )
                        median_secondary_entropy_vector.append(
                            np.median(secondary_entropy_data[i])
                        )
                    else:
                        mean_secondary_entropy_vector.append(-1)
                        median_secondary_entropy_vector.append(-1)

                # Format amplicon vectors
                amplicon_vectors = {
                    "mean_coverage_amplicon_vector": ", ".join(
                        [str(round(x, decimal_places)) for x in mean_coverage_vector]
                    )
                    if mean_coverage_vector
                    else "-",
                    "median_coverage_amplicon_vector": ", ".join(
                        [str(round(x, decimal_places)) for x in median_coverage_vector]
                    )
                    if median_coverage_vector
                    else "-",
                    "mean_entropy_amplicon_vector": ", ".join(
                        [str(round(x, decimal_places)) for x in mean_entropy_vector]
                    )
                    if mean_entropy_vector
                    else "-",
                    "median_entropy_amplicon_vector": ", ".join(
                        [str(round(x, decimal_places)) for x in median_entropy_vector]
                    )
                    if median_entropy_vector
                    else "-",
                    "mean_secondary_entropy_amplicon_vector": ", ".join(
                        [
                            str(round(x, decimal_places))
                            for x in mean_secondary_entropy_vector
                        ]
                    )
                    if mean_secondary_entropy_vector
                    else "-",
                    "median_secondary_entropy_amplicon_vector": ", ".join(
                        [
                            str(round(x, decimal_places))
                            for x in median_secondary_entropy_vector
                        ]
                    )
                    if median_secondary_entropy_vector
                    else "-",
                }

                # Display amplicon vectors
                for name, val in amplicon_vectors.items():
                    print(name, val, sep="\t")
