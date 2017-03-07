#!/usr/bin/env python

import sys
import os
import argparse
import numpy
from scipy.stats import rankdata

def parse_args(argv):
    ''' A method for taking in command line arguments and specifying
    help strings. '''
    parser = argparse.ArgumentParser(description="Parse a NetProphet adjacency matrix")
    parser.add_argument('-a', '--adjmtr', dest='adjmtr', type=str)
    parser.add_argument('-r', '--regulator', dest='regulator', type=str)
    parser.add_argument('-t', '--target', dest='target', type=str)
    parser.add_argument('-o', '--dir_output', dest='dir_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)
    parsed.dir_output = check_dir(parsed.dir_output)

    # load adjmtr, regulator and target lists
    adjmtr = numpy.loadtxt(parsed.adjmtr)
    tfs = numpy.loadtxt(parsed.regulator, dtype=str, delimiter="\t")
    targets = numpy.loadtxt(parsed.target, dtype=str, delimiter="\t")

    # write individual tf to target score files
    for i in range(adjmtr.shape[0]):
        writer = open(parsed.dir_output + tfs[i], "w")
        writer.write("#target\tscore\n")
        for j in range(adjmtr.shape[1]):
            if adjmtr[i, j] == 0:
                writer.write("%s\t0\n" % targets[j])
            else:
                writer.write("%s\t%0.17f\n" % (targets[j], adjmtr[i, j]))
        writer.close()

def check_dir(fd):
    if not fd.endswith('/'):
        fd += '/'
    if not os.path.exists(fd):
        os.makedirs(fd)
    return fd

if __name__ == "__main__":
    main(sys.argv)
