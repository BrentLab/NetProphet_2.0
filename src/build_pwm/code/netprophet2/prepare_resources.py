#!/usr/bin/python
import os
import sys
import argparse
import numpy as np

def parse_args(argv):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-g', '--genes', dest='genes')
    parser.add_argument('-r', '--regulators', dest='regulators')
    parser.add_argument('-e', '--expr_data', dest='expr_data')
    parser.add_argument('-f', '--fc_data', dest='fc_data')
    parser.add_argument('-c', '--conditions', dest='conditions')
    parser.add_argument('-or', '--output_reg_expr', dest='output_reg_expr')
    parser.add_argument('-of', '--output_data_fc', dest='output_data_fc')
    parser.add_argument('-oa', '--output_allowed', dest='output_allowed')
    parser.add_argument('-op1', '--output_pert_binary', dest='output_pert_binary')
    parser.add_argument('-op2', '--output_pert_boolean', dest='output_pert_boolean')
    # parser.add_argument('-ol', '--output_regulator_lists', dest='output_regulator_lists')
    parsed = parser.parse_args(argv[1:])
    return parsed


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


def write_tsv(fn, adjmtr, conds, genes):
	writer = open(fn, "w")
	for i in range(len(genes)):
		writer.write("%s\t" % genes[i])
	writer.write("\n")
	for i in range(len(adjmtr)):
		writer.write("%s\t" % conds[i])
		for j in range(len(adjmtr[i])):
			writer.write("%s\t" % adjmtr[i,j])
		writer.write("\n")
	writer.close()


def make_nonrepeat_conditions(conditions):
    conditions = np.array(conditions, dtype='|S500')
    non_unique_conds = [x for n, x in enumerate(conditions) if x in conditions[:n]]
    for x in np.unique(np.array(non_unique_conds)):
        indx = np.where(conditions == x)[0]
        for i in range(len(indx)):
            conditions[indx[i]] = str(conditions[indx[i]]) +'_'+ str(i+1)
    return conditions


def main(argv):
    parsed = parse_args(argv)
    ##make directories
    tmp_dir = os.path.dirname(os.path.abspath(parsed.output_reg_expr))
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    # if not os.path.exists(parsed.output_regulator_lists):
    #     os.makedirs(parsed.output_regulator_lists)

    ##load data
    genes = np.loadtxt(parsed.genes, dtype=str)
    regulators = np.loadtxt(parsed.regulators, dtype=str)
    conditions = np.loadtxt(parsed.conditions, dtype=str)

    ##prepare expression matrix of regulators (regulators x conditions)
    reg_indx = []
    for i in range(len(regulators)):
    	reg_indx.append(np.where(genes == regulators[i])[0][0])
    #data = np.loadtxt(parsed.expr_data)
    data = np.genfromtxt(parsed.expr_data, delimiter='\t') # changed by Dhoha.
    np.savetxt(parsed.output_reg_expr, data[reg_indx,:], fmt="%.10f", delimiter="\t")

    ##prepare fold change of expression data (conditions x genes)
    # data_fc = np.loadtxt(parsed.fc_data, dtype=str)
    data_fc = 2**data
    nonrepeat_conditions = make_nonrepeat_conditions(conditions)
    write_tsv(parsed.output_data_fc, data_fc.T, nonrepeat_conditions, genes)

    ##prepare allowed matrix (regulators x genes)
    allowed = np.ones((len(regulators),len(genes)), dtype=int)
    for i in range(len(reg_indx)):
    	allowed[i,reg_indx[i]] = 0
    np.savetxt(parsed.output_allowed, allowed, fmt="%d", delimiter="\t")

    ##prepare binary pert matrix (genes x conditions)
    ##and boolean pert matrix (conditions x genes)
    pert = np.zeros((len(genes),len(conditions)), dtype=int)
    pert_bool = np.array([['FALSE']*len(conditions)]*len(genes), dtype=str)
    for j in range(len(conditions)):
        if conditions[j] in regulators:
        	i = np.where(genes == conditions[j])[0]
        	if len(i) > 0:
        		pert[i[0],j] = 1
        		pert_bool[i[0],j] = 'TRUE'
    np.savetxt(parsed.output_pert_binary, pert, fmt="%d", delimiter="\t")
    write_tsv(parsed.output_pert_boolean, pert_bool.T, nonrepeat_conditions, genes)

    ##split regulators into sublists with 10 in each
    # step = 5
    # for i in range(int(np.ceil(len(regulators)/float(step)))):
    #     sublist = regulators[i*step:i*step+step]
    #     np.savetxt(parsed.output_regulator_lists+'/'+str(i+1)+'.txt', sublist, fmt="%s")

if __name__ == "__main__":
    main(sys.argv)
