from .api import all, query, iquery, parse_region
import argparse


def add_index(parser):
    parser.add_argument("--index", default=None)


def add_mapq(parser):
    parser.add_argument("--mapq", type=int, default=0)


def add_baseq(parser):
    parser.add_argument("--baseq", type=int, default=0)


def iterate(data, region=None):
    if region:
        chrom, start, end = parse_region(region)
        for (pos, ins_pos), row in sorted(data[chrom].items()):
            if (not start or pos >= start) and (not end or pos <= end):
                yield [pos, ins_pos, sum(row)] + row
    else:
        for chrom, chrom_data in data.items():
            for (pos, ins_pos), row in sorted(chrom_data.items()):
                yield [pos, ins_pos, sum(row)] + row


def run():
    parser = argparse.ArgumentParser()
    command = parser.add_subparsers(dest="command")

    all_parser = command.add_parser("all", help="Count bases at all positions.")
    all_parser.add_argument("bam")
    add_mapq(all_parser)
    add_baseq(all_parser)

    query_parser = command.add_parser(
        "query", help="Count bases in specific region (without an index file)."
    )
    query_parser.add_argument("bam")
    query_parser.add_argument("region")
    add_mapq(query_parser)
    add_baseq(query_parser)

    iquery_parser = command.add_parser(
        "iquery", help="Count bases in specific region (with an index file)."
    )
    iquery_parser.add_argument("bam")
    iquery_parser.add_argument("region")
    add_index(iquery_parser)
    add_mapq(iquery_parser)
    add_baseq(iquery_parser)

    args = parser.parse_args()

    print("\t".join(["pos", "ins", "cov", "a", "c", "g", "t", "ds", "n"]))

    if args.command == "all":
        data = all(args.bam, args.mapq, args.baseq)

        for row in iterate(data):
            print("\t".join(map(str, row)))

    elif args.command == "query":
        data = query(args.bam, args.region, args.mapq, args.baseq)

        for row in iterate(data, region=args.region):
            print("\t".join(map(str, row)))

    elif args.command == "iquery":
        data = iquery(args.bam, args.region, args.index, args.mapq, args.baseq)

        for row in iterate(data, region=args.region):
            print("\t".join(map(str, row)))

    else:
        print(
            """
██████╗░░█████╗░░██████╗███████╗███╗░░░███╗░█████╗░██████╗░
██╔══██╗██╔══██╗██╔════╝██╔════╝████╗░████║██╔══██╗██╔══██╗
██████╦╝███████║╚█████╗░█████╗░░██╔████╔██║███████║██████╔╝
██╔══██╗██╔══██║░╚═══██╗██╔══╝░░██║╚██╔╝██║██╔══██║██╔═══╝░
██████╦╝██║░░██║██████╔╝███████╗██║░╚═╝░██║██║░░██║██║░░░░░
╚═════╝░╚═╝░░╚═╝╚═════╝░╚══════╝╚═╝░░░░░╚═╝╚═╝░░╚═╝╚═╝░░░░░
        """
        )
