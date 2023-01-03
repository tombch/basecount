import pysam
import math
import argparse
import numpy as np


def open_samfile(bam):
    """
    Opens the file, temporarily suppressing HTSlib's messages while doing so
    """
    old_verbosity = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(bam, mode="rb")
    pysam.set_verbosity(old_verbosity)
    return samfile


def get_references(samfile, references=None):
    # Determine which references will be covered
    if references is None:
        # Calculate data on all references by default
        references = samfile.references
    else:
        references = set(references)

        # If references are provided, determine their validity
        for reference in references:
            if not (reference in samfile.references):
                raise Exception(f"{reference} is not a valid reference")

    # Put references into dictionary with their lengths as values
    return {ref: samfile.lengths[samfile.references.index(ref)] for ref in references}


def get_entropy(probabilities, normalised=False):
    entropy = sum([-(x * math.log2(x)) if x != 0 else 0 for x in probabilities])

    if normalised:
        return entropy / math.log2(len(probabilities))
    else:
        return entropy


def get_counts(samfile, references):
    base_counts = {}
    read_counts = {}

    for ref, ref_length in references.items():
        read_counts[ref] = 0

        for i in range(ref_length):
            base_counts[(ref, i, 0)] = {"A": 0, "C": 0, "G": 0, "T": 0, "DS": 0, "N": 0}

    for read in samfile.fetch(until_eof=True):
        # Reference name for read
        ref = read.reference_name

        if read.is_unmapped or (ref not in references):
            continue

        read_counts[ref] += 1

        # Index of where the read starts in the reference
        ref_pos = read.reference_start

        # Index for position in (the aligned portion of) the current read
        seq_pos = 0

        # The aligned portion of the read (soft clipped bases excluded)
        seq = read.query_alignment_sequence

        # Sequence of operations describing how the read is aligned to the reference
        # Each tuple is of the form (operation, length of operation)
        ctuples = read.cigartuples

        if (seq is not None) and (ctuples is not None):
            for operation, length in ctuples:
                # 0 : Matched (equal or not)
                # 7 : All bases equal to ref
                # 8 : All bases not equal to ref
                if operation == 0 or operation == 7 or operation == 8:
                    for r_p, s_p in zip(
                        range(ref_pos, ref_pos + length),
                        range(seq_pos, seq_pos + length),
                    ):
                        base_counts[(ref, r_p, 0)][seq[s_p]] += 1
                    seq_pos += length
                    ref_pos += length

                # 1 : Insertion in read
                elif operation == 1:
                    for i_p, s_p in zip(
                        range(1, 1 + length),
                        range(seq_pos, seq_pos + length),
                    ):
                        base_counts.setdefault(
                            (ref, ref_pos, i_p),
                            {"A": 0, "C": 0, "G": 0, "T": 0, "DS": 0, "N": 0},
                        )[seq[s_p]] += 1
                    seq_pos += length

                # 2 : Deletion in read
                # 3 : Skip in read
                elif operation == 2 or operation == 3:
                    for r_p in range(ref_pos, ref_pos + length):
                        base_counts[(ref, r_p, 0)]["DS"] += 1
                    ref_pos += length

                # Operations 4, 5 and 6 are soft clipping, hard clipping and padding respectively
                # None of these affect the read.query_alignment_sequence
                # Operation 9 is the 'back' operation which seems to be basically unheard of

    return base_counts, read_counts


class BaseCount:
    def __init__(
        self,
        bam,
        references=None,
        min_base_quality=0,
        min_mapping_quality=0,
        count_n_bases=False,
    ):
        self.columns = [
            "reference",
            "position",
            "insert_position",
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

        samfile = open_samfile(bam)
        self.references = get_references(samfile, references)
        self.base_counts, self.read_counts = get_counts(samfile, self.references)
        samfile.close()

        for ref, ref_length in self.references.items():
            for i in range(ref_length):
                count = list(self.base_counts[(ref, i, 0)].values())
                coverage = sum(count)
                probabilities = [
                    count / coverage if coverage > 0 else 0.0 for count in count
                ]
                percentages = [100 * probability for probability in probabilities]
                entropy = get_entropy(probabilities, normalised=True)
                secondary_count = list(count)
                secondary_count.pop(count.index(max(count)))
                secondary_coverage = sum(secondary_count)
                secondary_probabilities = [
                    count / secondary_coverage if secondary_coverage > 0 else 0.0
                    for count in secondary_count
                ]
                secondary_entropy = get_entropy(
                    secondary_probabilities, normalised=True
                )

                for k, v in zip(
                    [
                        "coverage",
                        "pc_a",
                        "pc_c",
                        "pc_g",
                        "pc_t",
                        "pc_ds",
                        "pc_n",
                        "entropy",
                        "secondary_entropy",
                    ],
                    [coverage] + percentages + [entropy, secondary_entropy],
                ):
                    self.base_counts[(ref, i, 0)][k] = v

                ins_pos = 1
                while self.base_counts.get((ref, i, ins_pos)):
                    count = list(self.base_counts[(ref, i, ins_pos)].values())
                    coverage = sum(count)
                    probabilities = [
                        count / coverage if coverage > 0 else 0 for count in count
                    ]
                    percentages = [100 * probability for probability in probabilities]
                    entropy = get_entropy(probabilities, normalised=True)
                    secondary_count = list(count)
                    secondary_count.pop(count.index(max(count)))
                    secondary_coverage = sum(secondary_count)
                    secondary_probabilities = [
                        count / secondary_coverage if secondary_coverage > 0 else 0
                        for count in secondary_count
                    ]
                    secondary_entropy = get_entropy(
                        secondary_probabilities, normalised=True
                    )

                    for k, v in zip(
                        [
                            "coverage",
                            "pc_a",
                            "pc_c",
                            "pc_g",
                            "pc_t",
                            "pc_ds",
                            "pc_n",
                            "entropy",
                            "secondary_entropy",
                        ],
                        [coverage] + percentages + [entropy, secondary_entropy],
                    ):
                        self.base_counts[(ref, i, ins_pos)][k] = v

                    ins_pos += 1

    def records(self, decimals=3):
        for ref, ref_length in self.references.items():
            for i in range(ref_length):
                record = self.base_counts[(ref, i, 0)]
                yield {
                    "reference": ref,
                    "position": i + 1,
                    "insert_position": 0,
                    "coverage": record["coverage"],
                    "num_a": record["A"],
                    "num_c": record["C"],
                    "num_g": record["G"],
                    "num_t": record["T"],
                    "num_ds": record["DS"],
                    "num_n": record["N"],
                    "pc_a": round(record["pc_a"], decimals),
                    "pc_c": round(record["pc_c"], decimals),
                    "pc_g": round(record["pc_g"], decimals),
                    "pc_t": round(record["pc_t"], decimals),
                    "pc_ds": round(record["pc_ds"], decimals),
                    "pc_n": round(record["pc_n"], decimals),
                    "entropy": round(record["entropy"], decimals),
                    "secondary_entropy": round(record["secondary_entropy"], decimals),
                }

                ins_pos = 1
                while self.base_counts.get((ref, i, ins_pos)):
                    record = self.base_counts[(ref, i, ins_pos)]
                    yield {
                        "reference": ref,
                        "position": i + 1,
                        "insert_position": ins_pos,
                        "coverage": record["coverage"],
                        "num_a": record["A"],
                        "num_c": record["C"],
                        "num_g": record["G"],
                        "num_t": record["T"],
                        "num_ds": record["DS"],
                        "num_n": record["N"],
                        "pc_a": round(record["pc_a"], decimals),
                        "pc_c": round(record["pc_c"], decimals),
                        "pc_g": round(record["pc_g"], decimals),
                        "pc_t": round(record["pc_t"], decimals),
                        "pc_ds": round(record["pc_ds"], decimals),
                        "pc_n": round(record["pc_n"], decimals),
                        "entropy": round(record["entropy"], decimals),
                        "secondary_entropy": round(
                            record["secondary_entropy"], decimals
                        ),
                    }
                    ins_pos += 1

    def mean_coverage(self, decimals=3):
        coverages = []
        for ref, ref_length in self.references.items():
            for i in range(ref_length):
                coverages.append(self.base_counts[(ref, i, 0)]["coverage"])
        return round(np.mean(coverages), decimals)

    def mean_entropy(self, decimals=3):
        entropies = []
        for ref, ref_length in self.references.items():
            for i in range(ref_length):
                entropies.append(self.base_counts[(ref, i, 0)]["entropy"])
        return round(np.mean(entropies), decimals)


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bam",
        help="Path to BAM file. An index file is not required",
    )
    parser.add_argument(
        "--refs",
        default=None,
        nargs="+",
        action="append",
        help="Choose specific reference(s) to run basecount on",
    )
    parser.add_argument(
        "--basequal",
        default=None,
        action="append",
        help="Minimum base quality. Default value: 0",
    )
    parser.add_argument(
        "--mapqual",
        default=None,
        action="append",
        help="Minimum mapping quality. Default value: 0",
    )
    parser.add_argument(
        "--decimals",
        default=None,
        action="append",
        help="How many decimal places to display. Default value: 3",
    )
    args = parser.parse_args()

    # Handle arguments
    references = handle_arg(args.refs, "refs")
    min_base_quality = int(handle_arg(args.basequal, "basequal", default=0, provided_once=True))  # type: ignore
    min_mapping_quality = int(handle_arg(args.mapqual, "mapqual", default=0, provided_once=True))  # type: ignore
    decimals = int(handle_arg(args.decimals, "decimals", default=3, provided_once=True))  # type: ignore

    bc = BaseCount(
        args.bam,
        references=references,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
    )

    print("\t".join(bc.columns))
    for x in bc.records(decimals=decimals):
        print("\t".join(map(str, x.values())))


if __name__ == "__main__":
    main()
