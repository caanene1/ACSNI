"""
ACSNI-derive
Author: Faraz Khan, PhD

"""
import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from ACSNI import dat, dbs
from ACSNI.corep import main_model_prediction


def cont_method(x, y):
    """
    Control methods {supervised}.

    Parameters
    ----------
    x: matrix
    y: target name

    Returns
    -------
    res: results matrix
    """

    xx = x.copy().T
    cols = dat.get_col_names(xx)
    xx["group"] = pd.cut(xx[y], bins=[np.min(xx[y]),
                                                   np.median(xx[y]),
                                                   np.max(xx[y])],
                            labels=["Low", "High"], include_lowest=True)

    res =  pd.DataFrame()
    for t in cols:
        x1 = np.array(xx[xx["group"] == "High"][t])
        y1 = np.array(xx[xx["group"] == "Low"][t])
        _, de = ttest_ind(x1, y1)

        de = np.array([de])
        frame1 = pd.DataFrame(de, index=[t])
        res = pd.concat([res, frame1], axis=0)
    res.rename(columns={0: "de.p"}, inplace=True)
    return res


def merge_minimal_derive(x):
    """
    Merge the minimal output results

    Parameters
    ----------
    x: Dictionary of results

    Returns
    -------
    end: Completion
    """

    mse = x[0]
    for i in x[1:]:
        mse = pd.merge(mse, i, how="outer", on="name")
    return mse

def get_ascni(prior_m, expression_m, mad, p, s, a, nn=".csv"):
    """
    Function to call ACSNI from derive.

    Parameters
    ----------
    prior_m: Prior matrix
    expression_m: Expression matrix
    mad: Mad cut
    p: percent
    s: seed
    a: p_value threshold for prediction call
    nn: The suffix to save output files

    Returns
    -------
    results

    """
    print("Searching")

    run_info = {
        "m": mad, "b": 1, "c": a,
        "p": p, "f": 0, "i": nn, "t": None,
        "w": None, "s": s}

    re = list()
    for i in prior_m:
        res = main_model_prediction(inp=expression_m,
                                        d_list=prior_m[i],
                                        gi=i, lp=p, s=s,
                                        run=run_info)

        re.append(res)

    re = merge_minimal_derive(re)
    return re

def get_cor(exp_m, cbtype, goi, madf, cort, corf, biotypefilter=False,
            biotypef=False,docor=False):
    """
    Function to correlate gene of interest against the genes of the same biotype

    Parameters
    ----------
    exp_m: Expression matrix
    biotypef: Gene Biotype file
    cbtype: Biotype of interest
    goi: Gene of interest
    madf: MAD threshold
    cort: Correlation threshold
    biotypefilter: If biotype file and filter parameters provided
    docor: Use modified correlation file provided
    corf: Modified correlation file to be used for the genration of prior 1.

    Returns
    -------
    gcor: Correlation matrix

    """

    if biotypefilter:
        if biotypef.columns[0] != 'gene' or biotypef.columns[1] != 'biotype':
            sys.exit("Make sure the column names and the order of the "
                     "biotype file is in the right format: 'gene','biotype'")
        else:
            bsub = biotypef[biotypef['biotype'].isin([cbtype])]
            esub = exp_m[exp_m['gene'].isin(bsub['gene'])]
    else:
        esub = exp_m

    if madf < 1:
        sys.exit("Choose MAD threshold >= 1")
    esub = dat.filter_uninformative(esub, madf)
    mtb=esub['gene'].str.contains(goi).sum()

    if mtb <= 0:
        sys.exit("Gene did not pass the MAD threshold. Try setting a "
        "lower value for -m/--mt, >= 1")
    esub = esub.dropna(how='any')
    esub.set_index('gene', inplace=True)
    esub.index.name = None

    if docor:
        gcor = corf
        gcor.sort_values(by=['cor'], inplace=True, ascending=False)
        gcor.index = gcor['gene']
        gcor.index.name = None
        gcor.to_csv("AC.csv", index=False, header=True)
        if cort < 0.6:
            sys.exit("Choose correlation threshold >= 0.6, non-negative")
        gcor = gcor[abs(gcor['cor']) > cort]

        if len(gcor.index) < 19:
            sys.exit("You don't have enough genes that show correlation "
            "above your set threshold. Try setting a lower correlation threshold")

        if len(gcor.index) > 50:
            gcor = gcor.iloc[0:50,:]
        gcor.to_csv("TopCC.csv", index=False, header=True)

        print("Finished extracting top correlated genes")
        return gcor

    else:
        gcor = esub.corrwith(esub.loc[goi], axis=1, drop=False, method='pearson')
        gde = cont_method(esub, goi)
        gcor = gcor.to_frame()
        gcor.columns = ['cor']
        gcor['gene'] = gcor.index
        gcor = gcor[['gene', 'cor']]
        cor_de = pd.concat([gcor, gde], axis=1)
        gcor.sort_values(by=['cor'], inplace=True, ascending=False)
        cor_de.sort_values(by=['cor'], inplace=True, ascending=False)
        cor_de.to_csv("AC.csv", index=False, header=True)

        if cort < 0.6:
            sys.exit("Choose correlation threshold >= 0.6")
        gcor = gcor[abs(gcor['cor']) > cort]
        if len(gcor.index) < 19:
            sys.exit("You don't have enough genes that show correlation "
                     "above your set threshold. Try a lower threshold")
        if len(gcor.index) > 50:
            gcor = gcor.iloc[0:50,:]
        gcor.to_csv("TopCC.csv", index=False, header=True)
        return gcor

def make_prior1(npc, gcor):
    """
    Function to make prior matrix for ACSNI-run 1

    Parameters
    ----------
    npc: Number of user specified prior matrix columns
    gcor: Correlation Matrix

    Returns
    -------
    prior1: Prior Matrix

    """

    cc = ['1', '2', '3', '4', '5']
    for i in range(6, npc + 1):
        cc.append(str(i))
    df = pd.DataFrame(np.random.randint(1, 2, size=(len(gcor.index),
                                                    int(len(cc)))), columns=cc)
    sd = 1
    ed = 3
    for i in cc:
        df.loc[sd:ed, i] = 0
        sd = ed + 1
        ed = ed + 3
    df1 = df.copy()
    df1.loc[0, :] = 0
    df.columns = ['V' + str(col) for col in df.columns]
    df1.columns = ['repV' + str(col) for col in df1.columns]
    prior1 = pd.concat([df, df1], axis=1)
    prior1.set_index(gcor.iloc[:, 0], inplace=True)
    prior1.index.name = None
    prior1.insert(0, 'gene', prior1.index)
    return prior1

def preprocess_run_one(pf, biotypefilter=False, biotypef=False, exclude=None):
    """
    Function to pre-process prediction output of ACSNI-run1

    Parameters
    ------------
    pf: Prediction output file from ACSNI-run1
    biotypef: Gene Biotype file
    biotypefilter: If biotype file and filter parameters provided.
    exclude: Biotype to exclude from Biotype file

    Returns
    -------
    predicf: Pre-processed prediction output

    """

    predicf = pd.DataFrame(data=pf)
    predicf.set_index('name', inplace=True)
    predicf.index.name = None
    predicf = predicf[predicf.sum(axis=1) >= 1]
    predicf.iloc[predicf > 0] = 1

    if biotypefilter:
        bsub2 = biotypef[~biotypef['biotype'].isin(exclude)]
        predicf = predicf[predicf.index.isin(bsub2['name'])]
    predicf.sort_index(inplace=True)
    return predicf

def process_run_one(x, y):
    """
    Function to find high confidence predicted genes from the pre-processed
    prediction output of ACSNI-run1.

    Parameters
    -----------
    x: ACSNI prediction results for option -f 0
    y: Length of initial prior

    Returns
    -------
    pdm: Processed results

    """

    pdm = list()

    for i in range(1, y + 1):
        m = "V" + str(i)
        r = "repV" + str(i)
        n = m + "vs" + r
        pdf = pd.DataFrame(x[m] > x[r], columns=[n])
        pdf = pdf[pdf[n].astype(str).str.contains("True")]
        pdm.append(pdf)

    mse = pdm[0]

    for i in pdm[1:]:
        mse = pd.concat([mse, i], axis=1)
    pdm = mse.fillna(0)
    pdm["sum"] = pdm.eq(True).sum(axis=1)
    pdm["sum"] = pdm["sum"].astype(int)
    pdm["prop"] = (pdm["sum"] / y) * 100
    pdm = pdm[pdm["prop"] >= 60]
    pdm = pdm.sort_values(by=["prop"], ascending=False)
    if len(pdm.index) > 80:
        pdm = pdm.iloc[0:80,:]
    return pdm

def make_prior2(hc, gt):
    """
    Function to make prior matrix 2 for ACSNI-run 2

    Parameters
    ----------
    hc: High-confidence predicted genes from ACSNI-run1
    gt: Gene of interest

    Returns
    -------
    prior2: Prior Matrix 2

    """
    cc2 = ['V1', 'V2']
    prior2 = pd.DataFrame(np.random.randint(1, 2, size=(len(hc.index),
    int(len(cc2)))), columns=cc2)
    prior2.set_index(hc.index, inplace=True)
    prior2.index.name = None
    goiindex = [gt]
    goi = pd.DataFrame([[1, 0]], columns=cc2, index=goiindex)
    prior2 = goi.append(prior2)
    prior2.insert(0, 'name', prior2.index)
    return prior2

def preprocess_run_two(pf2):
    """
    Function to pre-process prediction output of ACSNI-run2

    Parameters
    ------------
    pf2: Prediction output file from ACSNI-run2

    Returns
    -------
    predicf2: Pre-processed prediction output of ACSNI-run2

    """
    predicf2 = pd.DataFrame(data=pf2)
    predicf2.set_index('name', inplace=True)
    predicf2.index.name = None
    predicf2 = predicf2[predicf2.sum(axis=1) >= 1]
    predicf2.iloc[predicf2 > 0] = 1
    return predicf2

def process_run_two(hc2):
    """
    Function to find high confidence predicted genes from the pre-processed
    prediction output of ACSNI-run2.

    Parameters
    -----------
    hc2: High-confidence predicted genes from ACSNI-run2

    Returns
    -------
    pdm2: Processed results

    """
    pdm2 = pd.DataFrame(hc2['V1'] > hc2['V2'], columns=["V1vsV2"])
    pdm3 = pdm2
    pdm2 = pdm2[pdm2['V1vsV2'].astype(str).str.contains('True')]
    pdm2.to_csv("predictions.csv", index=True, header=True)
    pdm3.to_csv("FD.csv", index=True, header=True)
    return pdm2

def merge_multi(path):
    """
    Merge output from multiple De-novo genesets

    Parameters
    ----------
    path: Path to the working directory

    Returns
    -------
    out: List of dataframes D, N, co files
    """
    files = [f for f in os.listdir(path) if f.endswith(".csv")]
    nn, dd, co, ac, fd = list(), list(), list(), 0, 0

    for i in files:
        if i.startswith("N"):
            nn.append(pd.read_csv(os.path.join(path, i)))
            os.remove(os.path.join(path, i))
        elif i.startswith("D"):
            dd.append(pd.read_csv(os.path.join(path, i)))
            os.remove(os.path.join(path, i))
        elif i.startswith("AC"):
            ac = pd.read_csv(os.path.join(path, i))
            os.remove(os.path.join(path, i))
        elif i.startswith("FD"):
            fd = pd.read_csv(os.path.join(path, i))
            os.remove(os.path.join(path, i))
        elif i.startswith("code"):
            co.append(pd.read_csv(os.path.join(path, i)))
            os.remove(os.path.join(path, i))
        else:
           pass

    nn = [i.set_index('name') for i in nn]
    nn = pd.concat(nn, axis=1)

    dd = [i.set_index('name') for i in dd]
    dd = pd.concat(dd, axis=1)

    co = pd.concat(co, axis=1)

    return nn, dd, ac, fd, co

def save_merged_n_d_ac_fd_co(n, d, ac, fd, co, path, nfile, run_info):
    """
    Save multi-run

    Parameters
    ----------
    n: network
    d: direction
    ac: correlations
    fd: full result dataset
    co: code
    path: Path to saved file
    nfile: Name of the impute file
    run_info: Run arguments

    Returns
    -------
    """

    n.to_csv("N_{}".format(nfile), index=True)
    d.to_csv("D_{}".format(nfile), index=True)
    ac.to_csv("AC_{}".format(nfile), index=True)
    fd.to_csv("FD_{}".format(nfile), index=True)
    co.to_csv("code_{}".format(nfile), index=True)

    dbs_results = dbs.ACSNIDeriveResults(ac=pd.read_csv(os.path.join(path, "AC_{}".format(nfile))),
                                   n=pd.read_csv(os.path.join(path, "N_{}".format(nfile))),
                                   d=pd.read_csv(os.path.join(path, "D_{}".format(nfile))),
                                   fd=pd.read_csv(os.path.join(path, "FD_{}".format(nfile))),
                                   co=pd.read_csv(os.path.join(path, "code_{}".format(nfile))),
                                   run_info=run_info)

    for i in ["AC", "N", "D","FD","code"]:
        os.remove(os.path.join(path, "{}_{}".format(i, nfile)))

    mf = [f for f in os.listdir(path) if f.endswith(".h5")]
    for m in mf:
        os.remove(os.path.join(path, m))

    dbs.save_acsni(dbs_results, ("dbs" + nfile[:-4] + ".ptl"))
    return
