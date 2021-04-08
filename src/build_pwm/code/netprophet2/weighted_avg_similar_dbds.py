#!/usr/bin/python
import sys
import os
import glob
import argparse
import numpy as np

dbds_formats = ['multi_dbds', 'single_dbds']

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Merge TF network scores by weighted averaging, where the weight uses DBD-PWM similarity fit.")
    parser.add_argument("-n", "--fn_adjmtr_network", dest="fn_adjmtr_network")
    parser.add_argument("-r", "--fn_rids", dest="fn_rids")
    parser.add_argument("-a", "--dir_aligned_dbd", dest="dir_aligned_dbd")
    parser.add_argument("-d", "--dbd_cutoff", dest="dbd_cutoff", type=float, default=50)
    parser.add_argument("-f", "--dbds_formats", dest="dbds_formats", help="options: %s" % dbds_formats, default="single_dbds")
    parser.add_argument("-t", "--fn_dbd2rids_conversion", dest="fn_dbd2rids_conversion")
    parser.add_argument("-p", "--fn_pertrubed_rids", dest="fn_pertrubed_rids")
    parser.add_argument("-o", "--fn_output", dest="fn_output")
    parsed = parser.parse_args(argv[1:])
    return parsed


def sigmoid(x):
    return 0.9/(1+np.exp(-0.1*(x-40)))


def write_adjmtr(fn, adjmtr):
    writer = open(fn, "w")
    for i in range(len(adjmtr)):
        for j in range(len(adjmtr[i])):
            if adjmtr[i,j] == 0:
                writer.write("0\t")
            else:
                writer.write("%0.10f\t" % adjmtr[i,j])
        writer.write("\n")
    writer.close()


def get_regulators(fn_rids, fn_pert_rids):    
    rids = np.loadtxt(fn_rids, dtype=str)
    pert_rids = np.intersect1d(rids, np.loadtxt(fn_pert_rids, dtype=str)) if fn_pert_rids is not None else None
    return (rids, pert_rids)


def get_tf_weights_multi_dbds(dir_dbd, dbd_cutoff, fn_conv):
    ## parse tf-dbd conversion
    dbd2rid_list = np.loadtxt(fn_conv, dtype=str)
    dbd2rid_dict = {}
    for i in range(len(dbd2rid_list)):
        dbd2rid_dict[dbd2rid_list[i][0]] = dbd2rid_list[i][1]
    ## loop thru dbd similarity files
    tf_simscore_dict = {}
    for fn_dbd in glob.glob(dir_dbd + "*"):
        query_dbd = os.path.basename(fn_dbd)
        query_tf = dbd2rid_dict[query_dbd]
        if query_tf not in tf_simscore_dict.keys():
            tf_simscore_dict[query_tf] = {}
        ## get similarity scores 
        f = open(fn_dbd, "r")
        lines = f.readlines()
        f.close()
        for i in range(len(lines)):
            paired_tf = dbd2rid_dict[lines[i].split()[0]]
            paired_pctid = float(lines[i].split()[1])
            ## store all scores
            if paired_tf not in tf_simscore_dict[query_tf].keys():
                tf_simscore_dict[query_tf][paired_tf] = []
            tf_simscore_dict[query_tf][paired_tf].append(paired_pctid)
    ## use max tf similarity from multiple dbds to compute weights
    tf_weight_dict = {}
    for query_tf in tf_simscore_dict.keys():
        tf_weight_dict[query_tf] = {}
        for paired_tf, scores in tf_simscore_dict[query_tf].iteritems():
            scores = scores[scores >= dbd_cutoff]
            if len(scores) > 0:
                tf_weight_dict[query_tf][paired_tf] = sigmoid(max(scores))
    return tf_weight_dict


def get_tf_weights(dir_dbd, dbd_cutoff):
    ## loop thru dbd similarity files
    tf_simscore_dict = {}
    for fn_tf in glob.glob(dir_dbd + "*"):
        query_tf = os.path.basename(fn_tf)
        tf_simscore_dict[query_tf] = {}
        ## get similarity scores 
        f = open(fn_tf, "r")
        lines = f.readlines()
        f.close()
        for i in range(len(lines)):
            paired_tf = lines[i].split()[0]
            paired_pctid = float(lines[i].split()[1])
            ## store score
            tf_simscore_dict[query_tf][paired_tf] = paired_pctid
    ## use max tf similarity from multiple dbds to compute weights
    tf_weight_dict = {}
    for query_tf in tf_simscore_dict.keys():
        tf_weight_dict[query_tf] = {}
        for paired_tf, score in tf_simscore_dict[query_tf].iteritems():
            if score >= dbd_cutoff:
                tf_weight_dict[query_tf][paired_tf] = sigmoid(score)
    return tf_weight_dict


def update_tf_weights(tf_weight_dict, pert_rids):
    for query_tf in tf_weight_dict.keys():
        for paired_tf in tf_weight_dict[query_tf].keys():
            ## set non-pertrubed tf score to zero
            if (not paired_tf in pert_rids) and (paired_tf != query_tf):
                tf_weight_dict[query_tf][paired_tf] = 0


def average_scores(network_input, tf_weight_dict, rids):
    network_output = np.zeros(network_input.shape) 
    # weighted average
    for query_tf in rids:
        query_indx = np.where(rids == query_tf)[0][0]
        if (not query_tf in tf_weight_dict.keys()) or (len(tf_weight_dict[query_tf].keys()) < 2):
            ## unchanged score if query tf does not have other similar tfs
            network_output[query_indx, :] = network_input[query_indx, :]
        else:
            ## use allowed tfs 
            allowed_rids = np.intersect1d(tf_weight_dict[query_tf].keys(), rids)
            if len(allowed_rids) == 0:
                network_output[query_indx, :] = network_input[query_indx, :]
            else:
                paired_tf_indices = []
                paired_tf_weights = []
                for paired_tf in allowed_rids:
                    paired_tf_indices.append(np.where(rids == paired_tf)[0][0])
                    paired_tf_weights.append(tf_weight_dict[query_tf][paired_tf])
                # compute weighted average
                paired_tf_weights = np.array(paired_tf_weights)
                paired_tf_weights /= np.sum(paired_tf_weights)
                ## update the subnetwork
                network_sub = network_input[paired_tf_indices, :]
                network_output[query_indx, :] = np.dot(paired_tf_weights,  network_sub)
    return network_output


def main(argv):
    parsed = parse_args(argv)
    parsed.dir_aligned_dbd += "" if parsed.dir_aligned_dbd.endswith("/") else "/"

    ## get regualtor list
    rids, pert_rids = get_regulators(parsed.fn_rids, parsed.fn_pertrubed_rids)

    # parse dbd percent identity and compute weights from sigmoid fit
    sys.stdout.write("Computing TF weights ... ")
    if parsed.dbds_formats.lower() == "multi_dbds":
        if parsed.fn_dbd2rids_conversion is None:
            sys.exit("ERROR: DBD to TF conversion file must be provided if multi_dbds is used.")
        tf_weight_dict = get_tf_weights_multi_dbds(parsed.dir_aligned_dbd, parsed.dbd_cutoff, parsed.fn_dbd2rids_conversion)
    elif parsed.dbds_formats.lower() == 'single_dbds':
        tf_weight_dict = get_tf_weights(parsed.dir_aligned_dbd, parsed.dbd_cutoff)
    else:
        sys.exit("ERROR: Improper DBD alignment score format.")
    sys.stdout.write("DONE\n")

    # filter the weights for perturbed tfs only
    if pert_rids:
        sys.stdout.write("Recomputing TF weights for perturbed ones only ... ")
        update_tf_weights(tf_weight_dict, pert_rids)
        sys.stdout.write("DONE\n")

    # load network
    sys.stdout.write("Loading network ... ")
    network_input = np.loadtxt(parsed.fn_adjmtr_network)
    sys.stdout.write("DONE\n")

    # merge similar tfs using weighted average
    sys.stdout.write("Weighted averaging similar TFs ... ")
    network_output = average_scores(network_input, tf_weight_dict, rids)
    sys.stdout.write("DONE\n")

    # write weighted average network
    sys.stdout.write("Writing network ... ")
    write_adjmtr(parsed.fn_output, network_output)
    sys.stdout.write("DONE\n")


if __name__ == "__main__":
    main(sys.argv)
