#!/usr/bin/python
import sys
import argparse
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Map Network 1 to Network 2.")
    parser.add_argument('-n1', '--fn_net1', dest='fn_net1', type=str)
    parser.add_argument('-r1', '--fn_rids1', dest='fn_rids1', type=str)
    parser.add_argument('-g1', '--fn_gids1', dest='fn_gids1', type=str)
    parser.add_argument('-n2', '--fn_net2', dest='fn_net2', type=str)
    parser.add_argument('-r2', '--fn_rids2', dest='fn_rids2', type=str)
    parser.add_argument('-g2', '--fn_gids2', dest='fn_gids2', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)

    # read data
    adj1 = numpy.loadtxt(parsed.fn_net1, dtype=str)
    rids1 = numpy.loadtxt(parsed.fn_rids1, dtype=str)
    gids1 = numpy.loadtxt(parsed.fn_gids1, dtype=str)
    rids2 = numpy.loadtxt(parsed.fn_rids2, dtype=str)
    gids2 = numpy.loadtxt(parsed.fn_gids2, dtype=str)
    
    # map indices 
    rowinds = map_ind(rids1, rids2)
    colinds = map_ind(gids1, gids2)

    # build network 2
    writer = open(parsed.fn_net2, 'w')
    for i in range(len(rids2)):
    	for j in range(len(gids2)-1):
    		if rowinds[i] == -1 or colinds[j] == -1:
    			writer.write("0\t")
    		else:
    			score = adj1[rowinds[i],colinds[j]]
    			writer.write("%s\t" % score)
        if rowinds[i] == -1 or colinds[len(gids2)-1] == -1:
            writer.write("0\n")
        else:
            score = adj1[rowinds[i],colinds[len(gids2)-1]]
            writer.write("%s\n" % score)
    writer.close()

def map_ind(ids1, ids2):
    ind = [None] * len(ids2)
    for i in range(len(ids2)):
    	tempind = numpy.where(ids1 == ids2[i])[0]
    	if len(tempind) == 0:
    		ind[i] = -1
    	else:
    		ind[i] = tempind[0]
    return ind

if __name__ == "__main__":
    main(sys.argv)