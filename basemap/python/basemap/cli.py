from .api import all, query, iquery, parse_region
import argparse
import math


def add_index(parser):
    parser.add_argument("--index", default=None)


def add_mapq(parser):
    parser.add_argument("--mapq", type=int, default=0)


def add_baseq(parser):
    parser.add_argument("--baseq", type=int, default=0)


def add_stats(parser):
    parser.add_argument("--stats", action="store_true", default=False)


def get_entropy(probabilities, normalised=False):
    entropy = sum([-(x * math.log2(x)) if x != 0 else 0 for x in probabilities])

    if normalised:
        return entropy / math.log2(len(probabilities))
    else:
        return entropy


def get_stats(counts, dp=3):
    coverage = sum(counts)
    probabilities = [count / coverage if coverage > 0 else 0.0 for count in counts]
    percentages = [100 * probability for probability in probabilities]
    entropy = get_entropy(probabilities, normalised=True)
    secondary_count = list(counts)
    secondary_count.pop(counts.index(max(counts)))
    secondary_coverage = sum(secondary_count)
    secondary_probabilities = [
        count / secondary_coverage if secondary_coverage > 0 else 0.0
        for count in secondary_count
    ]
    secondary_entropy = get_entropy(secondary_probabilities, normalised=True)
    return (
        [coverage]
        + counts
        + [round(x, dp) for x in percentages + [entropy, secondary_entropy]]
    )


def iterate(data, region=None, stats=False):
    if region:
        chrom, start, end = parse_region(region)
        for (pos, ins_pos), row in sorted(data[chrom].items()):
            if (not start or pos >= start) and (not end or pos <= end):
                if stats:
                    yield [pos, ins_pos] + get_stats(row)
                else:
                    yield [pos, ins_pos, sum(row)] + row
    else:
        for chrom, chrom_data in data.items():
            for (pos, ins_pos), row in sorted(chrom_data.items()):
                if stats:
                    yield [pos, ins_pos] + get_stats(row)
                else:
                    yield [pos, ins_pos, sum(row)] + row


def run():
    parser = argparse.ArgumentParser()

    command = parser.add_subparsers(dest="command", required=True)

    all_parser = command.add_parser("all", help="Count bases at all positions.")
    all_parser.add_argument("bam")
    add_mapq(all_parser)
    add_baseq(all_parser)
    add_stats(all_parser)

    query_parser = command.add_parser(
        "query", help="Count bases in specific region (without an index file)."
    )
    query_parser.add_argument("bam")
    query_parser.add_argument("region")
    add_mapq(query_parser)
    add_baseq(query_parser)
    add_stats(query_parser)

    iquery_parser = command.add_parser(
        "iquery", help="Count bases in specific region (with an index file)."
    )
    iquery_parser.add_argument("bam")
    iquery_parser.add_argument("region")
    add_index(iquery_parser)
    add_mapq(iquery_parser)
    add_baseq(iquery_parser)
    add_stats(iquery_parser)

    args = parser.parse_args()

    if args.stats:
        print(
            "\t".join(
                [
                    "pos",
                    "ins",
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
                    "entropy",
                    "secondary_entropy",
                ]
            )
        )
    else:
        print("\t".join(["pos", "ins", "cov", "a", "c", "g", "t", "ds", "n"]))

    if args.command == "all":
        data = all(args.bam, args.mapq, args.baseq)

        for row in iterate(data, stats=args.stats):
            print("\t".join(map(str, row)))

    elif args.command == "query":
        data = query(args.bam, args.region, args.mapq, args.baseq)

        for row in iterate(data, region=args.region, stats=args.stats):
            print("\t".join(map(str, row)))

    elif args.command == "iquery":
        data = iquery(args.bam, args.region, args.index, args.mapq, args.baseq)

        for row in iterate(data, region=args.region, stats=args.stats):
            print("\t".join(map(str, row)))
