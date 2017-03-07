import sys
import os
import argparse
import numpy as nmp
import operator

"""
An importable file to provide utilities for model averaging in netprophet.
"""
def resort_by_weights(M, W):
    """ For all edges in M that have a corresponding edge in W,
    resort those edges in-place by their new score, M_ij*W_ij. 
    (or quantile normalization of multiplicated scores)"""
    # select only rows for which we weights to use
    W_avail = nmp.sum(W, 1) > 0
    working_net = nmp.abs(nmp.copy(M))
    # get the original values of the matrix where we have 
    orig_values = nmp.sort(nmp.ravel(nmp.copy(working_net)[W_avail]))[::-1]
    #working_net[W_avail] = (working_net[W_avail]+.001)*(W[W_avail]+.001)
    working_net[W_avail] = (working_net[W_avail])*(W[W_avail])
    index_to_value = []
    for j in range(nmp.shape(working_net)[0]):
        if not W_avail[j]:
            continue
        for i in range(nmp.shape(working_net)[1]):
            index_to_value.append([(j,i), working_net[j,i]])
    index_to_value.sort(key=lambda x: -x[1])
    for n,item in enumerate(index_to_value):
        j = item[0][0]; i = item[0][1]
        working_net[j, i] = orig_values[n]
        # working_net[j, i] = orig_values[n]*orig_signs[n]
    # working_net[W_avail] = working_net[W_avail]*orig_signs
    return working_net
    
def resort_by_pwm(network, pwm_net):
    """ """
    pwm_avail = nmp.sum(pwm_net, 1) > 0
    working_net = nmp.abs(nmp.copy(network))
    # get list of values from before pwms are used, map indices to values
    # orig_values = nmp.ravel(nmp.copy(working_net)[pwm_avail])
    # orig_signed = nmp.ravel(nmp.copy(network)[pwm_avail])
    # orig_signs[orig_signs!=0] = orig_signs[orig_signs!=0]/(nmp.abs(orig_signs[orig_signs!=0]))
    # orig_values = sorted(zip(list(orig_values), list(orig_signed)), key=lambda x:-x[0])
    # orig_signs = nmp.ravel(orig_signs)
    orig_values = nmp.sort(nmp.ravel(nmp.copy(working_net)[pwm_avail]))[::-1]
    working_net[pwm_avail] = (working_net[pwm_avail]+.001)*(pwm_net[pwm_avail]+.001)
    index_to_value = []
    for j in range(nmp.shape(working_net)[0]):
        if not pwm_avail[j]:
            continue
        for i in range(nmp.shape(working_net)[1]):
            index_to_value.append([(j,i), working_net[j,i]])
    index_to_value.sort(key=lambda x: -x[1])
    for n,item in enumerate(index_to_value):
        j = item[0][0]; i = item[0][1]
        working_net[j, i] = orig_values[n]
        # working_net[j, i] = orig_values[n]*orig_signs[n]
    # working_net[pwm_avail] = working_net[pwm_avail]*orig_signs
    return working_net
            
#end function

def list_geometric(ls):
    """ Returns the geometric mean of a list."""
    return reduce(operator.mul, ls)**(1.0/len(ls))

def rescale_matrix(mtr):
    ''' Given a matrix of float values, rescales it so the largest absolute
    value in the matrix is 1.'''
    return mtr/nmp.max(nmp.abs(mtr))
#end function

def rescale_shift_matrix(mtr):
    ''' Given a matrix of float values, shifts and rescales it so that the
    max absolute value of the matrix is 1, and the minimum absolute value
    of any value is 0 (i.e. the least negative is moved to 0, the most
    positive is 0).'''
    mtr_cp = nmp.copy(mtr)
    least_neg = nmp.max(mtr_cp[mtr_cp<0]) if nmp.any(mtr_cp<0) else 0
    mtr_cp[mtr_cp<0]-=least_neg
    least_pos = nmp.min(mtr_cp[mtr_cp>0]) if nmp.any(mtr_cp>0) else 0
    mtr_cp[mtr_cp>0]-=least_pos
    return rescale_matrix(mtr_cp)
#end function

def quadrant_combine(lasso_score, de_score, constants):
    ''' Returns a score for an edge based on its lasso score, de score, 
    and constants. '''
    cb_f = constants["Cb"]
    cd_f = constants["Cd"]
    w_f = None
    if lasso_score == 0 and de_score == 0:
        # if both scores are zero, we shouldn't give this edge any weight
        w_f = 0
    if lasso_score > 0 and de_score > 0:
        w_f = constants["quadrant I"]
    elif lasso_score < 0 and de_score > 0:
        w_f = constants["quadrant II"]
    elif lasso_score < 0 and de_score < 0:
        w_f = constants["quadrant III"]
    elif lasso_score > 0 and de_score < 0:
        w_f = constants["quadrant IV"]
    elif lasso_score != 0 and de_score == 0:
        w_f = constants["B"]
    elif lasso_score ==0 and de_score != 0:
        w_f = constants["D"]
    return (nmp.abs(lasso_score) + cb_f) * (nmp.abs(de_score) + cd_f) * w_f
#end function


def model_average_pwm_geometric(np_component, binding_strengths):
    ''' Performs geometric mean '''
    pwm_averaged = nmp.zeros(nmp.shape(np_component))
    for j, row in enumerate(np_component):
        for i, _ in enumerate(row):
            pwm_averaged[j, i] = list_geometric((np_component[j, i], 
                                                 binding_strengths[j, i]))
    pwm_averaged_nan_index = nmp.where(nmp.isnan(pwm_averaged))
    pwm_averaged[pwm_averaged_nan_index[0], pwm_averaged_nan_index[1]] = 0
    return pwm_averaged
#end function


def model_average_pwm_arithmetic(np_component, binding_strengths):
    ''' Performs arithmetic mean '''
    sign_np_component = nmp.sign(np_component)
    sign_np_component[nmp.where(sign_np_component == 0)] = 1
    abs_np_component = nmp.absolute(np_component)
    abs_binding_strengths = nmp.absolute(binding_strengths)
    combined = (abs_np_component + abs_binding_strengths)/2
    combined = nmp.multiply(combined, sign_np_component) 
    return combined
#end function


def model_average_pwm_arithmetic_intersect(np_component, binding_strengths):
    ''' Performs arithmetic mean of intersected edges '''
    sign_np_component = nmp.sign(np_component)
    sign_np_component[nmp.where(sign_np_component == 0)] = 1
    abs_np_component = nmp.absolute(np_component)
    abs_binding_strengths = nmp.absolute(binding_strengths)
    inds_intersect = nmp.nonzero(nmp.multiply(abs_np_component, abs_binding_strengths))
    combined = nmp.maximum(abs_np_component, abs_binding_strengths)
    combined[inds_intersect] = (abs_np_component[inds_intersect] + abs_binding_strengths[inds_intersect]) /2
    combined = nmp.multiply(combined, sign_np_component) 
    return combined
#end function


def model_average_np(lasso_component, de_component, 
                     constants = {"quadrant I":3, "quadrant II":1, 
                                  "quadrant III":1, "quadrant IV":1,
                                  "B":1, "D":2,
                                  "Cb":0.1, "Cd":0.01}):
    ''' Performs model averaging as in netprophet. '''
    betas = rescale_matrix(lasso_component)
    de = rescale_shift_matrix(de_component)
    retval = nmp.zeros(nmp.shape(lasso_component))
    for j, row in enumerate(lasso_component):
        if j%20 == 0:
            sys.stderr.write('working on row %i\n'%(j))
        for i, ix in enumerate(row):
            retval[j, i] = quadrant_combine(betas[j,i], 
                                            de[j,i],
                                            constants) 
    return retval
#end function


if __name__ == "__main__":
    raise NotImplemented("Module does not have CLI functionality.")
