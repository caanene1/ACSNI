"""
Name: ACSNI-data
Author: Chinedu A. Anene, Phd
"""

import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import numpy as np
import random, string, sys


def get_col_names(x):
    """
    Get names of numeric columns

    Parameters
    ----------
    x: Matrix
    """
    return x.select_dtypes([np.number]).columns

def get_row_index(x):
    """
    Get name of the row index.
    Throw error if there are more than categorical column.

    Parameters
    ----------
    x: Matrix
    """
    id_name = x.select_dtypes([np.object]).columns

    if len(id_name) != 1:
        sys.exit("ERROR: Make sure there just one ID column")

    return id_name

def parse_exp(x):
    """
    Process the expression file

    Parameters
    ----------
    x: matrix of expression

    Returns
    -------
    res: results matrix
    """

    x1 = get_row_index(x)
    x2 = get_col_names(x)
    y = pd.DataFrame()

    if len(x1) != 1:
        sys.exit("ERROR: Please, check the file and try again")

    y['name'] = x[x1[0]]
    y = pd.concat([y, x[x2]], axis=1)
    return y

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

    i_gene_len = len(df)
    if i_gene_len == 0:
        sys.exit("ERROR: Include column with gene name in the matrix")
    else:
        x = df.select_dtypes(include=["number"]).copy()
        row_deviance = x.sub(x.median(axis=1), axis=0)
        mad_row = np.abs(row_deviance).median(axis=1)
        temp = df[mad_row >= mad_f]

        if temp.shape[0] <= (i_gene_len/10):
            sys.exit(">90% uninformative genes. Please, "
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
    id_name = get_row_index(file_ss)
    scale = MinMaxScaler(feature_range=(0, 1))
    x = pd.DataFrame(scale.fit_transform(file_ss[cols]))
    out = pd.DataFrame(np.array(file_ss[id_name]))
    out.columns = ["name"]
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
    id_name = get_row_index(prior)
    prior = prior.rename(columns={id_name[0]: "name"})

    gene_set = {}
    for i in cols:
        gene_set[i] = pd.DataFrame(prior[prior[i] == 1]["name"])
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

    temp = pd.merge(get_scaled_values(input_file), d_list, how="inner", on="name")
    print("Geneset with {} genes in the expression".format(temp.shape[0]))

    if temp.shape[0] <= 11:
        sys.exit("Fewer than 12 genes in the set, use the f 2 option "
                 "to run existing method or check expression matrix")
    else:
        cols = get_col_names(temp)
        x = temp[cols]
        x = x.reindex(columns=sorted(x.columns))
        x = np.array(x)
        y = temp["name"]
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
    id_name = get_row_index(input_file)
    x = input_file[cols]
    x = x.reindex(columns=sorted(x.columns))
    y = input_file[id_name]
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

    nn = pd.DataFrame(y_train)
    w_nn = pd.merge(nn, w, how="left", on="name")
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