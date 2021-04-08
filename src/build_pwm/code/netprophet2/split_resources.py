#!/usr/bin/python
import sys
import argparse
import os.path
import numpy as np
from scipy.stats.mstats import gmean


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Split data on genes.")
    parser.add_argument('-i', '--dir_resources', type=str)
    parser.add_argument('-n', '--n_split', type=int, default=5)
    parsed = parser.parse_args(argv[1:])
    return parsed


def split_file(file_in, file_out_prefix, axis, n):
    ## axis: "c" - split columns; "r" - split rows
    ## special handling of tsv header
    is_tsv = file_in.endswith(".tsv")
    if is_tsv:
        data, rownames, colnames = load_tsv(file_in)
    else:
        data = np.loadtxt(file_in, dtype=str)
    is_vector = len(data.shape) == 1
    if is_vector:
        n_row = data.shape[0]
    else:
        n_row, n_col = data.shape

    ## split data
    step = int(n_row/float(n)) if axis == "r" else int(n_col/float(n))
    for i in range(n):
        if axis == "r":
            range_left = i*step
            range_right = max((i+1)*step, n_row) if i == n-1 else (i+1)*step
            data_out = data[range_left:range_right,]
            if is_tsv:
                rownames_out = rownames[range_left:range_right]
                colnames_out = colnames
        elif axis == "c":
            range_left = i*step
            range_right = max((i+1)*step, n_col) if i == n-1 else (i+1)*step
            data_out = data[:, range_left:range_right]
            if is_tsv:
                colnames_out = colnames[range_left:range_right]
                rownames_out = rownames
        ## save each split data
        file_out = "%s.%d" % (file_out_prefix, i+1)
        if is_tsv:
            write_tsv(file_out, (data_out, rownames_out, colnames_out))
        else:
            delim = "\n" if is_vector else "\t"
            np.savetxt(file_out, data_out, fmt="%s", delimiter=delim)


def load_tsv(file):
    cols = open(file, "r").readline().strip().split("\t")
    data = np.loadtxt(file, dtype=str, skiprows=1)
    rows = data[:,0]
    data = data[:,1:]
    return data, rows, cols


def write_tsv(file, data_tuple):
    data, rows, cols = data_tuple
    f = open(file, "w")
    f.write("%s\n" % "\t".join(cols))
    for i in range(len(rows)):
        f.write("%s\t%s\n" % (rows[i], "\t".join(data[i,:])))
    f.close()


def main(argv):
    parsed = parse_args(argv)
    # check arguments
    if not parsed.dir_resources.endswith("/"):
        parsed.dir_resources += "/"

    ## filename and split axis
    files_splits = [("data.expr", "r"),
                    ("genes", "r"),
                    ("de.signed.adj", "c"),
                    ("tmp/allowed.adj", "c"),
                    ("tmp/data.pert.adj", "r"),
                    ("tmp/data.fc.tsv", "c"),
                    ("tmp/data.pert.tsv", "c")]

    for filename, split_axis in files_splits:
        print "... working on %s" % filename
        file_in = parsed.dir_resources + filename
        file_out_prefix = parsed.dir_resources +"/"+ filename if filename.startswith("tmp") else parsed.dir_resources +"tmp/"+ filename
        split_file(file_in, file_out_prefix, split_axis, parsed.n_split)



if __name__ == "__main__":
    main(sys.argv)
