from pandas import DataFrame as Dff
import numpy as np
from ACSNI import dat
import pandas as pd
from scipy.stats import pearsonr, ttest_ind

"""
    elif fu == 2:
        corep.run_control(prior_m=prior_matrix[z], gi=z, 
        expression_m=expression_matrix, mad=mad_f,
                          nn=inp)
"""


def get_gs(xp, gs, mad):
    """
    Subset gene set expression for control methods

    Parameters
    ----------
    xp: expression matrix
    gs: Gene set
    mad: minimum absolute deviation

    Returns
    -------
    x: Expression of gene set
    """

    xp = dat.filter_uninformative(xp, mad_f=mad)
    temp = pd.merge(xp, gs, how="inner", on="gene")

    cols = dat.get_col_names(temp)
    x = temp[cols].T
    x = x.reindex(columns=sorted(x.columns))
    x = np.array(x)

    m_exp = Dff(np.median(x, axis=1))
    m_exp.rename(columns={0: "median"}, inplace=True)
    m_exp["group"] = pd.cut(m_exp["median"], bins=[np.min(m_exp["median"]),
                                                   np.median(m_exp["median"]),
                                                   np.max(m_exp["median"])],
                            labels=["Low", "High"], include_lowest=True)

    x_m = xp[dat.get_col_names(xp)].T
    x_m.rename(columns=xp["gene"], inplace=True)
    x_m.reset_index(drop=True, inplace=True)
    res = pd.concat([m_exp, x_m], axis=1)
    return res


def cont_methods(xp):
    """
    Control methods {unsupervised, supervised}

    Parameters
    ----------
    xp: expression matrix and derived

    Returns
    -------
    res: results matrix
    """

    cor_c = xp["median"]
    col = [col for col in xp.columns if col not in ["median", "group"]]

    cor_res = Dff()
    for i in col:
        x = xp[i]
        cr, _ = pearsonr(cor_c, x)
        cr = np.array([cr])
        frame = Dff(cr, index=[i])
        cor_res = pd.concat([cor_res, frame], axis=0)
    cor_res.rename(columns={0: "cor.r"}, inplace=True)

    de_res = Dff()
    for t in col:
        x = np.array(xp[xp["group"] == "High"][t])
        y = np.array(xp[xp["group"] == "Low"][t])
        _, de = ttest_ind(x, y)
        de = np.array([de])
        frame1 = Dff(de, index=[t])
        de_res = pd.concat([de_res, frame1], axis=0)
    de_res.rename(columns={0: "de.p"}, inplace=True)
    res = pd.merge(cor_res, de_res, left_index=True, right_index=True)
    return res