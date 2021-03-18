"""
Name: ACSNI-save
Author: Chinedu A. Anene, Phd
"""

from pandas import read_csv
import pandas as pd
import sys, os, random
from ACSNI.dbs import AcsniResults, save_acsni


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
        sys.exit("ERROR: Expression matrix must have numeric columns")

    new_exp = exp.select_dtypes(exclude=["number", "float", "int"]).copy()
    k_sam = len(new_exp.index)
    for i in nam_cols:
        new_exp[i] = random.sample(list(exp[i].values), k_sam)
    return new_exp

def save_merged_w_n_p_d(co, w, n, p, d, nn, run_info):
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
    run_info: Run parameters

    Returns
    -------
    """
    p_sum = p.select_dtypes(include=["number", "float", "int"]).copy()
    nnm = len(p_sum.columns)

    p["Sum_stat"] =  p_sum.astype(bool).sum(axis=1)
    p["Boot_Count"] = nnm

    dbs_results = AcsniResults(co=co, w=w, n=n, p=p, d=d,
                                   run_info=run_info)

    save_acsni(dbs_results, ("dbs" + nn[:-4] + ".ptl"))
    n1 = dbs_results.add_sub_process_direction(c=2)
    n1.to_csv("Network_{}".format(nn), index=False)
    return

def merge_multi(path, info, file):
    """
    Merge bootstrap run for summary start.
    Parameters
    ----------
    path: Path to the working directory
    info: Run information
    file: Input expression name

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
        elif i.startswith("W"):
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
        else :
            pass

    if len(co) == 1:
        co = co[0]

        ww = ww[0]
        ww.set_index('ID')

        nn = nn[0]
        nn.set_index('name')

        pp = pp[0]
        pp.set_index('name')

        dd = dd[0]
        dd.set_index('name')

    elif len(co) > 1:
        co = pd.concat(co, axis=0)

        ww = [i.set_index('ID') for i in ww]
        ww = pd.concat(ww, axis=1)

        nn = [i.set_index('name') for i in nn]
        nn = pd.concat(nn, axis=1)

        pp = [i.set_index('name') for i in pp]
        pp = pd.concat(pp, axis=1)

        dd = [i.set_index('name') for i in dd]
        dd = pd.concat(dd, axis=1)

    else:
        sys.exit("Fatal Error")

    save_merged_w_n_p_d(co=co, w=ww, n=nn, p=pp, d=dd, nn=file,
                        run_info=info)

    mf = [f for f in os.listdir(path) if f.endswith(".h5")]
    for m in mf:
        os.remove(os.path.join(path, m))

    return




