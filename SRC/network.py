"""
Name: ACSNI-network
Author: Chinedu A. Anene, Phd
"""

from SRC.dat import get_scaled_values
import pandas as pd
from sklearn.linear_model import LinearRegression,LogisticRegression
import numpy as np
import sys
from pandas import DataFrame as Dff

def unit_network(exp, w):
    """
    Construct unit wise regulatory networks

    Parameters
    ----------
    exp: Original expression matrix
    w: Model generated sample weights

    Returns
    -------
    lm_w: Matrix of the association between functional unit and predicted genes
    """

    w1 = w.select_dtypes(include=["number", "float", "int"]).copy()
    exp1 = get_scaled_values(exp)

    exp_set = exp1.select_dtypes(include=["number", "float", "int"]).copy().T

    lm_w = exp1[["gene"]]
    lm_w1 = exp1[["gene"]]
    for i in range(w1.shape[1]):
        x = np.array(w1.iloc[:, i]).reshape((-1, 1))

        lm_g = []
        lm_d = []
        for t in range(exp_set.shape[1]):
            y = np.array(exp_set.iloc[:, t]).reshape((-1, 1))
            lm_r = LinearRegression(n_jobs=-1)
            lm_r.fit(x, y)
            lm_g.append(lm_r.score(x, y))
            lm_d.append(float(lm_r.coef_))
        lm_g = pd.DataFrame(lm_g)
        lm_d = pd.DataFrame(lm_d)
        lm_g.columns.values[0] = i
        lm_d.columns.values[0] = i
        lm_w = pd.concat([lm_w, lm_g], axis=1)
        lm_w1 = pd.concat([lm_w1, lm_d], axis=1)
    return lm_w, lm_w1

def beta_params(mu, var):
    """
    Estimate the parameters of beta distribution

    Parameters
    ----------
    mu: Mean of x variable
    var: variance of x variable

    Returns
    -------
    estimates: the beta and alpha
    """
    alpha = ((1- mu) / var - 1 / mu) * mu ** 2
    beta = alpha * (1 / mu - 1)
    return alpha, beta

def get_explained_variation(exp, nm, cla_ss):
    """
    Get level of interaction between phenotype and subprocesses.

    Parameters
    ----------
    exp: dbs object from get_associated_module
    nm: phenotype variable name
    cla_ss: Type of variable

    Returns
    -------
    Frame: The variation level

    """
    weights = exp.select_dtypes(include=["number", "float", "int"]).columns
    lm_w = Dff({"sub":weights})

    lm_g = []
    for i in weights:
        x = np.array(exp[i]).reshape(-1, 1)
        y = exp[[nm]].values.ravel()

        if cla_ss == "character":
            lm_r = LogisticRegression(n_jobs=-1)
        elif cla_ss == "numeric":
            lm_r = LinearRegression(n_jobs=-1)
        else:
            sys.exit("Invalid value set for parameter -c. Use numeric or character for the type of phenotype")

        lm_r.fit(x, y)
        lm_r.score(x, y)
        lm_g.append(lm_r.score(x, y))
    lm_g = Dff(lm_g)
    lm_w = pd.concat([lm_w, lm_g], axis=1)
    return lm_w

def get_summary_stat_phenotype(x, nm):
    """
    Phenotype association statistics
    Parameters
    ----------
    x: Subprocess vs phenotype association table
    nm: phenotype variable name

    Returns
    -------
    x: Update table and print

    """
    qq = x[[0]].quantile([.25, .75])
    q25 = qq.loc[0.25, 0]
    q75 = qq.loc[0.75, 0]
    #
    m = float(x[[0]].mean())
    sd = float(x[[0]].std())
    print("Statistics of variations in subprocesses explained by {}".format(nm))
    print("q25 {} \n q75 {} \n mean {} \n std {}".format(q25, q75, m, sd))
    x["Association"] = np.where(x[0] >= q75, "Strong", "Weak")
    return x
