"""
Name: ACSNI-core
Author: Chinedu A. Anene, Phd
"""

import sys
import pandas as pd
import numpy as np
from pandas import DataFrame as Dff
from scipy.stats import beta
from SRC import mud, dat
from SRC.network import unit_network, beta_params



def check_input(nc):
    """
    Check number of samples

    Parameters
    ----------
    nc: Number of columns

    Returns
    -------
    end: Pass or early stop
    """

    print("Checking number of samples")
    if nc <= 15:
        sys.exit("{} samples is small for the tool"
                 .format(nc))
    elif nc >= 600:
        sys.exit("{} samples is large, use "
                 "ACSNI-split".format(nc))
    else:
        print("{} samples is suitable".format(nc))
    return

def check_name(e_df):
    """
    Check expression matrix

    Parameters
    ----------
    e_df: Expression matrix

    Returns
    -------
    end: Pass or early exit
    """
    c_col = len(e_df.columns)
    num_exp = e_df.select_dtypes(include=["number"]).copy()
    n_col = c_col - len(num_exp.columns)

    if n_col >= 2:
        print("{} categorical columns found".format(n_col))
        sys.exit("Remove the {} extra label columns in "
                 "the expression matrix".format(n_col-1))
    else:
        print("Completed expression matrix checks with 0 exit")
    return

def check_duplicates(df):
    """
    Check duplicated genes

    Parameters
    ----------
    df: Expression matrix

    Returns
    -------
    end: Pass or early exit
    """
    is_dup_row = df.duplicated(subset="gene", keep='first')
    dup_genes = df[is_dup_row]["gene"]

    if dup_genes.size >= 1:
        print("Found {} duplicated genes".format(dup_genes.size))
        print(dup_genes)
        sys.exit("Remove duplicated genes and try again")
    else:
        print("Completed duplicate checks with 0 exit")
    return

def score_class(x, p_cut):
    """
    Call interactions {"P", "B} and counts them for global prediction

    Parameters
    ----------
    x: LM Scores and genes
    p_cut: probability threshold

    Returns
    -------
    x: Updated score frame with "Predicted" columns"
    """

    print("Predicting context-specific network")
    nn = dat.get_col_names(x)

    for i in nn:
        a, b = beta_params(x[i].mean(), x[i].var())

        x[i] = 1 - beta.cdf(x[i], a, b)
        x[i] = np.where(x[i] <= p_cut, "P", "B")
    x["Predicted"] = x.iloc[:].eq("P").sum(axis=1).to_frame()
    return x

def main_model_prediction(inp, w, d_list, mad, full_r, i, nn, lp, s, p_cut):
    """
    Run final model

    Parameters
    ----------
    inp: expression
    w: Optional weights
    d_list: Prior
    mad: minimum absolute deviation
    full_r: mode
    i: Gene set name
    s: Boot number
    nn: Input name
    lp: proportion for layer
    p_cut: p_value threshold for calling gene significant

    Returns
    -------
    count_pre: predictions
    weights: functional units
    """
    predict, weight = (0,0)

    exp_mat = dat.filter_uninformative(inp, mad_f=mad)
    x_train, y_train = dat.scale_first(input_file=exp_mat, d_list=d_list)
    check_input(x_train.shape[1])

    if w is None:
        w_nn = w
    else:
        w_nn = dat.get_model_weights(w, y_train)

    weight = mud.mud_output(x_train, w_nn, lp)
    weight = Dff(weight)
    weight["ID"] = inp.select_dtypes(include=["number", "float", "int"]).columns
    un, ud = unit_network(exp_mat, weight)

    print("Saving network (N) for {} code".format(i))
    un.to_csv("N{}_{}_{}".format(s, i, nn), index=False)

    print("Saving network (D) for {} code".format(i))
    ud.to_csv("D{}_{}_{}".format(s, i, nn), index=False)

    if full_r == 1:
        predict = score_class(un, p_cut)

    elif full_r == 0:
        predict = score_class(un, p_cut)
        predict = predict[[col for col in predict.columns if col in ["gene", "Predicted"]]]
    else:
        print("Please set the correct argument -f, one of 1, 0, 2")
    return predict, weight

def run_full(prior_m, gi, expression_m, w, mad, f, nn, p, s, a):
    """
    Full mode and outputs

    Parameters
    ----------
    prior_m: Gene set matrix
    gi: Gene set name
    expression_m: Full expression matrix
    w: Optional weight matrix
    mad: Minimum median absolute deviation
    f: Mode
    nn: Name of input expression matrix
    p: Percentage for layer size
    s: Boot number
    a: P_value threshold for calling gene

    Returns
    -------
    end: outputs all results to the current working directory
    """
    count_pre, weights = main_model_prediction(inp=expression_m, w=w,
                                               d_list=prior_m,
                                               mad=mad, full_r=f,
                                               i=gi, nn=nn, lp=p, s=str(s), p_cut=a)


    print("Saving sample weights (W) for {}-{}".format(gi, str(s)))
    weights.to_csv("W{}_{}.csv".format(str(s), gi), index=False)

    print("Saving prediction (P) for {}-{} code".format(gi,str(s)))
    count_pre.to_csv("P{}_{}_{}".format(str(s), gi, nn), index=False)
    return

def run_minimal(prior_m, gi, expression_m, w, mad, f, p, s, a, nn=0):
    """
    Minimal output mode

    Parameters
    ----------
    prior_m: Prior matrix
    gi: Gene set name
    expression_m: Expression matrix
    w: Optional weights
    mad: Median absolute deviation
    f: Mode
    p: Percentage for layer size
    s: Boot number
    a: P_value threshold for calling gene

    Returns
    -------
    re: Dictionary of Gene set results
    """
    count_pre, _ = main_model_prediction(inp=expression_m, w=w, d_list=prior_m,
                                         mad=mad, full_r=f, i=gi, nn=nn, lp=p, s=s, p_cut=a)

    count_pre.columns = ["gene", gi]
    return count_pre

def merge_minimal_out(x, nn):
    """
    Merge the minimal output results

    Parameters
    ----------
    x: Dictionary of results
    nn: Name of input file

    Returns
    -------
    end: Completion
    """

    mse = x[0]
    for i in x[1:]:
        mse = pd.merge(mse, i, how="outer", on="gene")
    mse.to_csv("predicted_{}.csv".format(nn), index=False)
    return
