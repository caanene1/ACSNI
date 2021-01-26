"""
Name: ACSNI-data
Author: Chinedu A. Anene, Phd
"""

import pandas as pd
from pandas import DataFrame as Dff
from sklearn.preprocessing import MinMaxScaler
import numpy as np
import random, string, sys


def get_col_names(file_s):
    """
    Get numeric columns

    Parameters
    ----------
    file_s: Expression matrix with gene names

    Returns
    -------
    names: Names of numeric columns
    """

    return [col for col in file_s.columns if col not in ["gene"]]

def remove_unexpressed(df):
    """
    Remove unexpressed genes

    Parameters
    ----------
    df: Full expression matrix

    Returns
    -------
    res: Filtered expression matrix
    """

    x = df.select_dtypes(include=["number", "float", "int"]).copy()
    nc = x.shape[1]
    x['count'] = x.T[x.T == 0].count()
    res = df[x['count'] <= nc]
    res.reset_index(drop=True, inplace=True)
    return res

def filter_uninformative(df, mad_f):
    """
    Filter uninformative genes

    Parameters
    ----------
    df: Expression matrix
    mad_f: Minimum median absolute deviation

    Returns
    -------
    temp: Filtered expression matrix
    """

    print("Filtering uninformative genes by MAD")
    i_gene_len = len(df)
    if i_gene_len == 0:
        sys.exit("Please, include column with gene name in the matrix")
    else:
        x = df.select_dtypes(include=["number"]).copy()
        row_deviance = x.sub(x.median(axis=1), axis=0)
        mad_row = np.abs(row_deviance).median(axis=1)

        temp = df[mad_row >= mad_f]
        if temp.shape[0] <= (i_gene_len/4):
            sys.exit(">25% uninformative genes. Please, "
                     "check the expression matrix or reduce -m. "
                     "Use pre-filtered matrix if appropriate")
        else:
            temp.reset_index(drop=True, inplace=True)
    return temp

def get_scaled_values(file_ss):
    """
    Scale the expression on samples

    Parameters
    ----------
    file_ss: Filtered expression matrix

    Returns
    -------
    res: Scaled expression matrix
    """

    cols = get_col_names(file_ss)
    scale = MinMaxScaler(feature_range=(0, 1))
    x = Dff(scale.fit_transform(file_ss[cols]))
    out = Dff(np.array(file_ss["gene"]))
    out.columns = ["gene"]
    res = pd.concat([out, x], axis=1)
    return res

def gene_sets(prior):
    """
    Process the prior matrix

    Parameters
    ----------
    prior: Gene set matrix of 1s and 0s

    Returns
    -------
    gene_set: Dictionary of prior matrix
    """

    cols = get_col_names(prior)
    gene_set = {}
    for i in cols:
        gene_set[i] = pd.DataFrame(prior[prior[i] == 1]["gene"])
    return gene_set

def scale_first(input_file, d_list):
    """
    Subset gene set expression

    Parameters
    ----------
    input_file: Filtered, scaled expression matrix
    d_list: Single gene set names

    Returns
    -------
    x: Gene set expression
    y: Gene names
    """

    print("Processing the prior expression profiles")
    temp = pd.merge(get_scaled_values(input_file), d_list, how="inner", on="gene")
    print("Found {} genes from the set in the expression matrix".format(temp.shape[0]))

    if temp.shape[0] <= 11:
        sys.exit("Fewer than 12 genes in the set, use the f 2 option "
                 "to run existing method or check expression matrix")
    else:
        cols = get_col_names(temp)
        x = temp[cols]
        x = x.reindex(columns=sorted(x.columns))
        x = np.array(x)
        y = temp["gene"]
    return x, y

def scale_second(input_file):
    """
    Process expression for prediction

    Parameters
    ----------
    input_file: Complete expression matrix

    Returns
    -------
    y: Gene names
    x: Expression values
    """

    cols = get_col_names(input_file)
    x = input_file[cols]
    x = x.reindex(columns=sorted(x.columns))
    y = input_file["gene"]
    return y, x

def get_model_weights(w, y_train):
    """
    Prepare optional weight matrix

    Parameters
    ----------
    w: weight matrix
    y_train: Training gene names (Gene set)

    Returns
    -------
    w_nn1: Formatted weight matrix
    """

    nn = Dff(y_train)
    w_nn = pd.merge(nn, w, how="left", on="gene")
    w_nn["w"] = w_nn["w"].fillna(0)
    w_nn1 = w_nn["w"]
    return w_nn1

def name_generator(l=6, s=string.ascii_uppercase + string.digits):
    """
    Function to generate random names for saving files.
    :param l: Length of the string.
    :param s: Where to sample from.
    :return:
    """
    return ''.join(random.choice(s) for _ in range(l))