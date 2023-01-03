import argparse
import basemap


def add_index(parser):
    parser.add_argument("--index", default=None)


def add_mapq(parser):
    parser.add_argument("--mapq", type=int, default=0)


def add_baseq(parser):
    parser.add_argument("--baseq", type=int, default=0)


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

    if args.command == "all":
        data = basemap.all(args.bam, args.mapq, args.baseq)

        pos = 1
        for chrom in data.values():
            while chrom.get((pos, 0)):
                print("\t".join(map(str, [pos] + chrom[(pos, 0)])))
                pos += 1

    elif args.command == "query":
        region, start, end = basemap.parse_region(args.region)

        data = basemap.query(args.bam, args.mapq, args.baseq, args.region)

        if start:
            pos = start
        else:
            pos = 1

        while data[region].get((pos, 0)) and not (end and pos > end):
            print("\t".join(map(str, [pos] + data[region][(pos, 0)])))
            pos += 1

    elif args.command == "iquery":
        if not args.index:
            args.index = args.bam + ".bai"

        region, start, end = basemap.parse_region(args.region)

        data = basemap.iquery(
            args.bam, args.bam + ".bai", args.mapq, args.baseq, args.region
        )

        if start:
            pos = start
        else:
            pos = 1

        while data[region].get((pos, 0)) and not (end and pos > end):
            print("\t".join(map(str, [pos] + data[region][(pos, 0)])))
            pos += 1

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
