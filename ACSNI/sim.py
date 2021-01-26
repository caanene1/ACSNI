"""
Name: ACSNI-save
Author: Chinedu A. Anene, Phd
"""

import random
from pandas import read_csv
import pandas as pd
import sys
import os
from ACSNI import dbs


def shuffle_exp(exp):
    """
    Shuffle expression matrix

    Parameters
    ----------
    exp: Original expression matrix

    Returns
    -------
    lm_w: Shuffled expression matrix
    """

    nam_cols = exp.select_dtypes(include=["number", "float", "int"]).columns
    if len(nam_cols) == 0:
        sys.exit("Expression matrix must have numeric columns")

    new_exp = exp.select_dtypes(exclude=["number", "float", "int"]).copy()
    k_sam = len(new_exp.index)
    for i in nam_cols:
        new_exp[i] = random.sample(list(exp[i].values), k_sam)
    return new_exp

def merge_multi(path):
    """
    Merge bootstrap run for summary start
    Parameters
    ----------
    path: Path to the working directory

    Returns
    -------
    out: List of dataframes W, N, P, Co
    """
    files = [f for f in os.listdir(path) if f.endswith(".csv")]
    ww, nn, pp, dd, co = list(), list(), list(), list(), list()

    for i in files:
        if i.startswith("code"):
            co.append(read_csv(os.path.join(path, i)))
            os.remove(os.path.join(path, i))
        if i.startswith("W"):
            ww.append(read_csv(os.path.join(path, i)))
            os.remove(os.path.join(path, i))
        elif i.startswith("N"):
            nn.append(read_csv(os.path.join(path, i)))
            os.remove(os.path.join(path, i))
        elif i.startswith("P"):
            pp.append(read_csv(os.path.join(path, i)))
            os.remove(os.path.join(path, i))
        elif i.startswith("D"):
            dd.append(read_csv(os.path.join(path, i)))
            os.remove(os.path.join(path, i))
        elif i.endswith("_model.h5"):
            os.remove(os.path.join(path, i))
        else:
            pass

    if len(co) > 1:
        co = pd.concat(co, axis=0)

        ww = [i.set_index('ID') for i in ww]
        ww = pd.concat(ww, axis=1)

        nn = [i.set_index('gene') for i in nn]
        nn = pd.concat(nn, axis=1)

        pp = [i.set_index('gene') for i in pp]
        pp = pd.concat(pp, axis=1)

        dd = [i.set_index('gene') for i in dd]
        dd = pd.concat(dd, axis=1)
    return co, ww, nn, pp, dd

def save_merged_w_n_p_d(co, w, n, p, d, nn, path, run_info):
    """
    Save multi-run
    Parameters
    ----------
    w: weights
    n: networks
    p: predictions
    d: direction
    co: Pathway code
    nn: Name of input expression matrix
    path: Path to saved file
    run_info: Run parameters

    Returns
    -------
    """
    p_sum = p.select_dtypes(include=["number", "float", "int"]).copy()
    nnm = len(p_sum.columns)

    p["Sum_stat"] =  p_sum.astype(bool).sum(axis=1)
    p["Boot_Count"] = nnm

    co.to_csv("code_{}".format(nn), index=True)
    w.to_csv("W_{}".format(nn), index=True)
    n.to_csv("N_{}".format(nn), index=True)
    d.to_csv("D_{}".format(nn), index=True)
    p.to_csv("P_{}".format(nn), index=True)

    dbs_results = dbs.AcsniResults(co=read_csv(os.path.join(path, "code_{}".format(nn))),
                                   w=read_csv(os.path.join(path, "W_{}".format(nn))),
                                   n=read_csv(os.path.join(path, "N_{}".format(nn))),
                                   p=read_csv(os.path.join(path, "P_{}".format(nn))),
                                   d=read_csv(os.path.join(path, "D_{}".format(nn))),
                                   run_info=run_info)

    for i in ["code", "D", "W"]:
        os.remove(os.path.join(path, "{}_{}".format(i, nn)))


    dbs.save_acsni(dbs_results, ("dbs" + nn[:-4] + ".ptl"))
    p1 = dbs_results.get_global_prediction()
    n1 = dbs_results.add_sub_process_direction(2)

    n1.to_csv("N_{}".format(nn), index=False)
    p1.to_csv("P_{}".format(nn), index=False)
    return


