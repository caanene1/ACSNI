"""
Name: ACSNI-core
Author: Chinedu A. Anene, Phd
"""

import sys
import pandas as pd
from ACSNI.mud import DimRed
from ACSNI.dat import get_row_index, filter_uninformative, \
    scale_first, get_model_weights
from ACSNI.network import NetScore


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

    if nc <= 15:
        sys.exit("{} samples is small for the tool"
                 .format(nc))
    elif nc >= 800:
        sys.exit("{} samples is large, use "
                 "ACSNI-split".format(nc))
    else:
        print("{} samples".format(nc))
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
        print("ERROR: {} categorical columns found".format(n_col))
        sys.exit("Remove the {} extra label columns in "
                 "the expression matrix".format(n_col-1))
    else:
        pass
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

    id_name = get_row_index(df)

    is_dup_row = df.duplicated(subset=id_name, keep='first')
    dup_genes = df[is_dup_row][id_name]

    if dup_genes.size >= 1:
        print("ERROR: Found {} duplicated genes".format(dup_genes.size))
        print(dup_genes)
        sys.exit("Remove duplicated genes and try again")
    else:
        pass
    return


def main_model_prediction(inp, d_list, gi, lp, s, run):
    """
    Run final model

    Parameters
    ----------
    inp: expression
    d_list: Prior
    gi: Gene set name
    s: Boot number
    lp: Layer dimension
    run: The details of the run to be performed

    Returns
    -------
    count_pre: predictions
    weights: functional units
    aede: optimal dimension
    """

    exp_mat = filter_uninformative(inp, mad_f=run.get("m"))
    x_train, y_train = scale_first(input_file=exp_mat, d_list=d_list)
    check_input(x_train.shape[1])

    if run.get("w") is None:
        w_nn = run.get("w")

    else:
        w_nn = get_model_weights(pd.read_csv(run.get("w")), y_train)

    mod = DimRed(x=x_train, w = w_nn, p = lp)
    mod.fit()
    mod.add_reduced_row(inp.select_dtypes(include=["number","float","int"]).columns)

    net = NetScore(x=exp_mat, w=mod.get_reduced(), p=run.get("c"))
    net.fit()
    net.save_net(str(s), gi, run.get("i"))

    if run.get("f") == 1:
        mod.get_reduced().to_csv("W{}_{}.csv".format(str(s), gi), index=False)

        net.get_predicted(run.get("f")).to_csv("P{}_{}_{}".format(str(s), gi,
                                                                  run.get("i")), index=False)
        count_p = mod.get_aede()

    elif run.get("f") == 0:
        count_p = net.get_predicted(run.get("f"))
        count_p.columns = ["name", gi]

    else:
        sys.exit("ERROR: Invalid value set for parameter -f")

    return count_p

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
        mse = pd.merge(mse, i, how="outer", on="name")
    mse.to_csv("predicted_{}.csv".format(nn), index=False)
    return
