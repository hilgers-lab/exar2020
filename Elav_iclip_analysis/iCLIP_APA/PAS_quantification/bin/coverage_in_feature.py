#
# @author: Barbara Hummel
# @modified by Michael Rauer
#

#!/usr/bin/env python
import argparse
import pyBigWig
import numpy as np


def sum_coverage(bigwig_file, bed_file):
    counts = {}
    bigwig = pyBigWig.open(bigwig_file)

    for line in open(bed_file):
        if line.startswith('#'):
            continue

        cols = line.strip().split()
        # print(line)
        # strand specific signal
        vals = bigwig.values(cols[0], int(cols[1]), int(cols[2]))
        countSum = np.nansum(vals)
        # save in counts
        if cols[3] in counts:
            c = counts[cols[3]][4]
            counts[cols[3]] = [cols[0], cols[1], cols[2], cols[3], c + countSum, cols[5]]
        else:
            counts[cols[3]] = [cols[0], cols[1], cols[2], cols[3], countSum, cols[5]]
    bigwig.close()

    return counts

def parseArgs():
    parser = argparse.ArgumentParser(description="Counts number of 5' ends in each region of a BED-file.")
    parser.add_argument("bigwig", help="Bigwig file")
    parser.add_argument("bedfile", help="Input bed file")
    parser.add_argument("output", help="Output file coverage read counts for each region.")
    return parser

def main():
    args = parseArgs().parse_args()

    readCounts = sum_coverage(args.bigwig, args.bedfile)

    o = open(args.output, "w")
    o.write("Chromosome\tStart\tEnd\tGene\tCounts\tStrand\n")
    for region in readCounts:
        o.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            readCounts[region][0],
            readCounts[region][1],
            readCounts[region][2],
            readCounts[region][3],
            int(readCounts[region][4]),
            readCounts[region][5]))
    o.close()


if __name__ == "__main__":
    main()




