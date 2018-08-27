#!/usr/bin/python
import sys
import argparse
import os.path
import numpy as np
from scipy.stats.mstats import gmean


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Build motif network using the (inferred) PWM scores on promoters.")
    parser.add_argument('-i', '--fn_inferred', dest='fn_inferred', type=str)
    parser.add_argument('-r', '--fn_rids', dest='fn_rids', type=str)
    parser.add_argument('-g', '--fn_gids', dest='fn_gids', type=str)
    parser.add_argument('-f', '--dir_fimo', dest='dir_fimo', type=str)    
    parser.add_argument('-o', '--fn_adjmtr', dest='fn_adjmtr', type=str)
    parser.add_argument('-s', '--summary_suffix', dest='summary_suffix', type=str, default=".summary")
    parser.add_argument('-t', '--thld_type', dest='thld_type', type=str)
    parser.add_argument('-v', '--thld_val', dest='thld_val', type=float, default=0)
    parsed = parser.parse_args(argv[1:])
    return parsed


def build_network(fn_rids, fn_gids, fn_inferred, dir_fimo, summary_suffix, thld_type, thld_val):
    # initialize network
    rids = np.loadtxt(fn_rids, dtype=str)
    gids = np.loadtxt(fn_gids, dtype=str)
    adjmtr = np.zeros([len(rids), len(gids)])
    # iterate thru motifs
    f = open(fn_inferred, "r")
    lines = f.readlines()
    f.close()
    for line in lines:
        # get inferred tf and database motif
        inferred, _, _, zscore, robust = line.strip().split('\t')
        motifs = inferred.split(',')

        if thld_type is None:
            ## all inferred motifs are valid
            indx = np.where(rids == inferred)[0]
            if len(indx) > 0:
                continue
            ## get fimo score and build subnetwork
            subadjmtr = build_subnetwork(motifs, gids, dir_fimo, summary_suffix)
            adjmtr[indx[0], :] = gmean(subadjmtr).data

        else:
            # get confidence score and threshold type
            zscore = float(zscore)
            robust = float(robust.split("/")[0])
            if thld_type == 'zscore':
                score = zscore
            elif thld_type == 'robust':
                score = robust
            ## check if inferred motif passes motif quality threshold
            if score < thld_val:
                continue
            indx = np.where(rids == inferred)[0]
            if len(indx) <= 0:
                continue
            ## get fimo score and build subnetwork
            subadjmtr = build_subnetwork(motifs, gids, dir_fimo,summary_suffix)
            if len(motifs) > 1:
                adjmtr[indx[0], :] = gmean(subadjmtr).data
            else:
                adjmtr[indx[0], :] = subadjmtr
    return adjmtr


def build_subnetwork(motifs, gids, dir_fimo, summary_suffix):
    adjmtr = np.zeros([len(motifs), len(gids)])
    for j in range(len(motifs)):
        fn_motif = dir_fimo + motifs[j] + summary_suffix
        if not os.path.isfile(fn_motif):
            continue
        score_dict = get_fimo_scores(fn_motif)
        for k in range(len(gids)):
            gid = gids[k]
            adjmtr[j, k] = score_dict[gid] if gid in score_dict else 0
    return adjmtr


def get_fimo_scores(fn):
    ## load file
    f = open(fn, "r")
    lines = f.readlines()
    f.close()
    ## store target scores into dictionary
    d = {}
    for line in lines:
        line = line.strip().split()
        name = line[1]
        score = max([float(line[3]), float(line[5])])
        d[name] = score
    return d


def write_adjmtr(adjmtr, fn):
    writer = open(fn, "w")
    for i in range(len(adjmtr)):
        for j in range(len(adjmtr[i])):
            if adjmtr[i,j] == 0:
                writer.write("0.\t")
            else:
                writer.write("%0.10f\t" % adjmtr[i,j])
        writer.write("\n")
    writer.close()


def main(argv):
    parsed = parse_args(argv)
    # check arguments
    if not parsed.dir_fimo.endswith("/"):
        parsed.dir_fimo += "/"
    if parsed.thld_type is not None:
        if parsed.thld_type not in ["zscore", "robust"]:
            sys.exit("No confidence threshold to filter motifs.\n")

    ## build network
    sys.stdout.write("Building motif network ... ")
    network = build_network(parsed.fn_rids, parsed.fn_gids, parsed.fn_inferred, parsed.dir_fimo, parsed.summary_suffix, parsed.thld_type, parsed.thld_val)
    # write adjmtr file
    sys.stdout.write("DONE\nWriting network ... ")
    write_adjmtr(network, parsed.fn_adjmtr)
    sys.stdout.write("DONE\n")


if __name__ == "__main__":
    main(sys.argv)
