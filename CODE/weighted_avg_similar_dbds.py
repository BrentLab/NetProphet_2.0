#!/usr/bin/python
import sys
import os
import glob
import argparse
import numpy

dbds_formats = ['multi_dbds', 'single_dbds', 'single_adjmtr']

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

def main(argv):
    parsed = parse_args(argv)
    parsed.dir_aligned_dbd += "" if parsed.dir_aligned_dbd.endswith("/") else "/"

    rids = numpy.loadtxt(parsed.fn_rids, dtype=str)

    use_pert_rids_only = False
    if parsed.fn_pertrubed_rids != None:
        use_pert_rids_only = True
        pert_rids = numpy.intersect1d(rids, numpy.unique(numpy.loadtxt(parsed.fn_pertrubed_rids, dtype=str)))

    # parse dbd percent identity and compute weights from sigmoid fit
    sys.stdout.write("computing TF weights ...\n")

    if parsed.dbds_formats.lower() == "multi_dbds":
        # parse tf-dbd conversion
        dbd2rid_list = numpy.loadtxt(parsed.fn_dbd2rids_conversion, dtype=str)
        dbd2rid_dict = {}
        for i in range(len(dbd2rid_list)):
            dbd2rid_dict[dbd2rid_list[i][0]] = dbd2rid_list[i][1]

        # get max tf similarity from multiple dbds and compute weights
        tf_fractions_dict = {}
        fns_dbd = glob.glob(parsed.dir_aligned_dbd + "*")
        for fn_dbd in fns_dbd:
            query_dbd = os.path.basename(fn_dbd)
            query_tf = dbd2rid_dict[query_dbd]

            f = open(fn_dbd, "r")
            lines = f.readlines()
            for i in range(len(lines)):
                paired_tf = dbd2rid_dict[lines[i].split()[0]]
                paired_pctid = float(lines[i].split()[1])

                if paired_pctid >= parsed.dbd_cutoff:
                    paired_weight = sigmoid(paired_pctid)

                    if query_tf in tf_fractions_dict.keys():
                        if (not paired_tf in tf_fractions_dict[query_tf].keys()) or (paired_pctid > tf_fractions_dict[query_tf][paired_tf]):
                            tf_fractions_dict[query_tf][paired_tf] = paired_weight
                    else:
                        tf_fractions_dict[query_tf] = {paired_tf: paired_weight}
            f.close()

    elif parsed.dbds_formats.lower() == 'single_dbds':
        # parse single dbd similarity and compute tf weigths
        tf_fractions_dict = {}
        fns_tf = glob.glob(parsed.dir_aligned_dbd + "*")
        for fn_tf in fns_tf:
            query_tf = os.path.basename(fn_tf)

            f = open(fn_tf, "r")
            lines = f.readlines()
            for i in range(len(lines)):
                paired_tf = lines[i].split(":")[0]
                paired_pctid = float(lines[i].split()[1])

                if paired_pctid >= parsed.dbd_cutoff:
                    paired_weight = sigmoid(paired_pctid)

                    if query_tf in tf_fractions_dict.keys():
                        if (not paired_tf in tf_fractions_dict[query_tf].keys()) or (paired_pctid > tf_fractions_dict[query_tf][paired_tf]):
                            tf_fractions_dict[query_tf][paired_tf] = paired_weight
                    else:
                        tf_fractions_dict[query_tf] = {paired_tf: paired_weight} 
            f.close()

    elif parsed.dbds_formats.lower() == 'single_adjmtr':
        pass

    else:
        sys.exit("Improper dbd alignment score format.")

    # filter the weights for perturbed tfs only
    if use_pert_rids_only:
        sys.stdout.write("recomputing TF weights for perturbed ones only\n")

        for query_tf in tf_fractions_dict.keys():
            for paired_tf in tf_fractions_dict[query_tf].keys():
                if (not paired_tf in pert_rids) and (paired_tf != query_tf):
                    tf_fractions_dict[query_tf][paired_tf] = 0

    # merge similar tfs using weighted average
    sys.stdout.write("weighted averaging similar TFs ...\n")

    # parse network
    network_input = numpy.loadtxt(parsed.fn_adjmtr_network)
    network_output = numpy.zeros(network_input.shape) 

    # weighted average
    for query_tf in rids:
        query_tf_index = numpy.where(rids == query_tf)[0][0]

        if (not query_tf in tf_fractions_dict.keys()):
            network_output[query_tf_index, :] = network_input[query_tf_index, :]

        else:
            if (len(tf_fractions_dict[query_tf].keys()) > 1):
                paired_tf_indices = numpy.array([], dtype=int)
                paired_tf_weights = numpy.array([], dtype=int)
                for paired_tf in numpy.intersect1d(tf_fractions_dict[query_tf].keys(), rids):
                    paired_tf_indices = numpy.append(paired_tf_indices, numpy.where(rids == paired_tf)[0][0])
                    paired_tf_weights = numpy.append(paired_tf_weights, tf_fractions_dict[query_tf][paired_tf])

                # compute weighted average
                paired_tf_weights /= numpy.sum(paired_tf_weights)
                network_sub = network_input[paired_tf_indices, :]
                network_output[query_tf_index, :] = numpy.dot(paired_tf_weights,  network_sub)
            else:
                network_output[query_tf_index, :] = network_input[query_tf_index, :]

    # write weighted average network
    sys.stdout.write("writing network ...\n")

    write_adjmtr(parsed.fn_output, network_output)


def sigmoid(x):
    return 0.9/(1+numpy.exp(-0.1*(x-40)))


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


if __name__ == "__main__":
    main(sys.argv)
