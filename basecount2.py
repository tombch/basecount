from basecount.main import open_samfile, get_references
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    args = parser.parse_args()

    samfile = open_samfile(args.bam)
    references = get_references(samfile)
    references = {
        ref: samfile.lengths[samfile.references.index(ref)] for ref in references
    }

    counts = {}

    for ref, ref_length in references.items():
        for i in range(ref_length):
            counts[(ref, i, 0)] = {"A": 0, "C": 0, "G": 0, "T": 0, "DS": 0, "N": 0}

    for read in samfile.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        # Reference name for read
        ref = read.reference_name

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
                        counts[(ref, r_p, 0)][seq[s_p]] += 1
                    seq_pos += length
                    ref_pos += length

                # 1 : Insertion in read
                elif operation == 1:
                    for i_p, s_p in zip(
                        range(1, 1 + length),
                        range(seq_pos, seq_pos + length),
                    ):
                        counts.setdefault(
                            (ref, ref_pos, i_p),
                            {"A": 0, "C": 0, "G": 0, "T": 0, "DS": 0, "N": 0},
                        )[seq[s_p]] += 1
                    seq_pos += length

                # 2 : Deletion in read
                # 3 : Skip in read
                elif operation == 2 or operation == 3:
                    for r_p in range(ref_pos, ref_pos + length):
                        counts[(ref, r_p, 0)]["DS"] += 1
                    ref_pos += length

                # Operations 4, 5 and 6 are soft clipping, hard clipping and padding respectively
                # None of these affect the read.query_alignment_sequence
                # Operation 9 is the 'back' operation which seems to be basically unheard of

    print(
        "\t".join(
            [
                "ref",
                "pos",
                "ins_pos",
                "cov",
                "a",
                "c",
                "g",
                "t",
                "ds",
                "n",
                "pc_a",
                "pc_c",
                "pc_g",
                "pc_t",
                "pc_ds",
                "pc_n",
            ]
        )
    )
    for ref, ref_length in references.items():
        for i in range(ref_length):
            coverage = sum(counts[(ref, i, 0)].values())
            print(
                "\t".join(
                    map(
                        str,
                        [ref, i + 1, 0, coverage]
                        + list(counts[(ref, i, 0)].values())
                        + [
                            round(100 * x / coverage, 3) if coverage > 0 else 0
                            for x in counts[(ref, i, 0)].values()
                        ],
                    )
                )
            )

            ins_pos = 1
            while counts.get((ref, i, ins_pos)):
                coverage = sum(counts[(ref, i, 0)].values())
                print(
                    "\t".join(
                        map(
                            str,
                            [ref, i + 1, ins_pos, coverage]
                            + list(counts[(ref, i, ins_pos)].values())
                            + [
                                round(100 * x / coverage, 3) if coverage > 0 else 0
                                for x in counts[(ref, i, ins_pos)].values()
                            ],
                        )
                    )
                )
                ins_pos += 1


if __name__ == "__main__":
    main()
