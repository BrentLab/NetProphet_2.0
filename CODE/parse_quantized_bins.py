#!/usr/bin/python

import sys
import os
import argparse
import operator
import glob
import os.path

def parse_args(argv):
    parser = argparse.ArgumentParser(description=
        "Discretize continuous NetProphet scores into specific number of equal-value-range bins, along with a bin of all zeros.")
    parser.add_argument('-i', '--dir_input', dest='dir_input', type=str)
    parser.add_argument('-n', '--num_bin', dest='num_bin', type=int, default=20)
    parser.add_argument('-o', '--dir_output', dest='dir_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed


def errprint(st):
    sys.stderr.write(st + "\n")


def check_dir(fd):
    if not fd.endswith('/'):
        fd += '/'
    return fd


def process_score(fi, fo, num_bin):
    # parse scores and put (target: score) in dictionary
    dict_target = {}
    lines = open(fi).readlines()
    for i, line in enumerate(lines):
        if i > 0 and line.strip():
            line_split = line.split();
            dict_target[line_split[0]] = float(line_split[1])
    # sort dictionary based on score
    dict_target_sorted = sorted(dict_target.items(), key=operator.itemgetter(1))
    dict_target_sorted.reverse()

    # get score intervals
    #print(dict_target_sorted)
    try:
        max_score = dict_target_sorted[0][1]
    except IndexError:
        print('dict_target_sorted: ', dict_target_sorted)
        exit()
              
    min_score = dict_target_sorted[len(dict_target_sorted)-1][1]
    interval_score = (max_score-min_score)/(num_bin-1)
    
    list_target_sorted = []
    # do binning on non-zero scores
    i = 0
    j = 0
    while i < len(dict_target_sorted):
        if round(dict_target_sorted[i][1],15) >= round((max_score-interval_score*(j+1)),15):
            if dict_target_sorted[i][1] != 0:
                list_target_sorted.append((dict_target_sorted[i][0], j, dict_target_sorted[i][1]))
            i += 1
        else:   j += 1
    # do binning on zero scores
    for i in range(len(dict_target_sorted)):
        if dict_target_sorted[i][1] == 0:
            list_target_sorted.append((dict_target_sorted[i][0], num_bin-1, dict_target_sorted[i][1]))

    # write discrete tf-target scores to output file
    writer = open(fo, "w")
    writer.write("target\tscore\n")
    for i in range(len(list_target_sorted)):
        writer.write("%s\t%d\n" % (list_target_sorted[i][0], list_target_sorted[i][1]))
    writer.close()


def main(argv):
    parsed = parse_args(argv)
    parsed.dir_input = check_dir(parsed.dir_input)
    parsed.dir_output = check_dir(parsed.dir_output)
    if not os.path.exists(parsed.dir_output):
        os.makedirs(parsed.dir_output)

    fns = glob.glob(parsed.dir_input + "*")
    for fn in fns:
        tf = os.path.basename(fn)
        process_score(fn, parsed.dir_output + tf, parsed.num_bin)


if __name__ == "__main__":
    main(sys.argv)
