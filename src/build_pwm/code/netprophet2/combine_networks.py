import sys
import os
import argparse
from model_averaging_utils import *
import numpy as nmp


"""
Code to take the LASSO and DE components of netprophet and aggregate them.

This module will make use of PWM information to suss out what TF-target
interactions are likely or not.
"""

averaging_strategies = {'NP':model_average_np, 
                        'geometric':model_average_pwm_geometric,
                        'arithmetic':model_average_pwm_arithmetic,
                        'arithmetic_intersect':model_average_pwm_arithmetic_intersect,
                        'resort':resort_by_weights}

def parse_args(argv):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-l','--lasso_component', dest='lasso_component')
    parser.add_argument('-d','--de_component', dest='de_component')
    parser.add_argument('-n','--np_component', dest='np_component')
    parser.add_argument('-od','--output_dir', dest='output_dir')
    parser.add_argument('-om','--output_adjmtr_name', dest='output_adjmtr_name')
    parser.add_argument('-ol','--output_adjlst_name', dest='output_adjlst_name')
    parser.add_argument('-r','--regulator_names', dest='regulator_names')
    parser.add_argument('-b','--binding_strengths', dest='binding_strengths',
                        default = None)
    parser.add_argument('-g','--target_names', dest='target_names')
    parser.add_argument('-s','--strategy', dest='strategy', default='NP',
                        help='options: %s'%str(averaging_strategies.keys()))
    parsed = parser.parse_args(argv[1:])
    return parsed

def output(result):
    """ Write the results out the stdout. """
    pass

def main(argv):
    """ The main module should take in inputs from the command line, 
    carry out the 'CLI' functionality of the script, then write the
    results to stdout. """
    parsed = parse_args(argv)
    
    assert parsed.strategy in averaging_strategies, "provided strategy must " + \
        "be one of the following: %s"%(str(averaging_strategies))
    
    # check input network arguments
    if (parsed.lasso_component == None and parsed.de_component == None and parsed.np_component != None):
        is_np_net = True
    elif (parsed.lasso_component != None and parsed.de_component != None and parsed.np_component == None):
        is_np_net = False
    else:
        sys.exit("Input network argument error.")

    # read in output directory
    if (not is_np_net):

        sys.stderr.write("Reading input arguments ... ")

        # read in LASSO values
        lasso_component = nmp.loadtxt(parsed.lasso_component)
        # read in DE values
        de_component = nmp.loadtxt(parsed.de_component)
        # read in desired name for final combined adjmtr
        output_adjmtr_name = os.path.join(parsed.output_dir, 
                                          parsed.output_adjmtr_name)
        # read in PWM binding information, if available
        binding_strengths = None
        if parsed.binding_strengths != None:
            binding_strengths = nmp.loadtxt(parsed.binding_strengths)
            
        # Optional:
        making_adjlst = False
        if (parsed.output_adjlst_name != None and 
            parsed.regulator_names != None and
            parsed.target_names != None):
            
            making_adjlst = True
            # read in desired name for final combined adjlst
            output_adjlst_name = os.path.join(parsed.output_dir,
                                              parsed.output_adjlst_name)
            # # read in list of regulator names
            # regulator_names = [row.strip() for row in open(parsed.regulator_names)]
            # # read in list of target gene names
            # target_names = [row.strip() for row in open(parsed.target_names)]

        sys.stderr.write("Done\nCombining with motif network ... ")
        
        # perform model averaging
        if parsed.strategy == 'NP':
            combined = model_average_np(lasso_component, de_component)
        else:
            np_component = model_average_np(lasso_component, de_component)
            combined = averaging_strategies[parsed.strategy](np_component, binding_strengths)

        sys.stderr.write("Done\n")

    else:

        sys.stderr.write("Reading input arguments ... ")

        # read in np values
        np_component = nmp.loadtxt(parsed.np_component)
        # read in desired name for final combined adjmtr
        output_adjmtr_name = os.path.join(parsed.output_dir, 
                                          parsed.output_adjmtr_name)
        # read in PWM binding information, if available
        binding_strengths = None
        if parsed.binding_strengths != None:
            binding_strengths = nmp.loadtxt(parsed.binding_strengths)
            
        # Optional:
        making_adjlst = False
        if (parsed.output_adjlst_name != None and 
            parsed.regulator_names != None and
            parsed.target_names != None):
            
            making_adjlst = True
            # read in desired name for final combined adjlst
            output_adjlst_name = os.path.join(parsed.output_dir,
                                              parsed.output_adjlst_name)

        sys.stderr.write("Done\nCombining with motif network ... ")
        
        # perform model averaging
        if parsed.strategy == 'NP':
            combined = np_component
        else:
            combined = averaging_strategies[parsed.strategy](np_component,binding_strengths)

        sys.stderr.write("Done\n")
    
    # if args provided, write out combined lists as an adjacency list
    sys.stderr.write("Writing adjmtr file ... ")
    # nmp.savetxt(output_adjmtr_name, combined)
    write_adjmtr(output_adjmtr_name, combined)
    sys.stderr.write("Done\n")
#end function

def write_adjmtr(fn, adjmtr):
    writer = open(fn, "w")
    for i in range(len(adjmtr)):
        for j in range(len(adjmtr[i])):
            if adjmtr[i,j] == 0:
                writer.write("0\t")
            else:
                writer.write("%0.15f\t" % adjmtr[i,j])
        writer.write("\n")
    writer.close()

if __name__ == "__main__":
    main(sys.argv)
