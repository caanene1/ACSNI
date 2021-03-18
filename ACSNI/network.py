"""
Name: ACSNI-network
Author: Chinedu A. Anene, Phd
"""

from ACSNI.dat import get_scaled_values
import pandas as pd
from sklearn.linear_model import LinearRegression
from scipy.stats import beta
import numpy as np
from ACSNI.dat import get_col_names


class NetScore:
    """
    Class for network construction and scoring.
    """

    def __init__(self, x, w, p):
        """
        Network Class.

        Parameters:
            x: Input matrix
            w: Weights to merge and calculate
            p: probability threshold
        """
        self.x = x
        self.w = w
        self.p = p
        self.__un = None
        self.__ud = None

    def __str__(self):
        return "Network construction class"

    def __beta_params(self, x):
        """
        Estimate the parameters of beta distribution.

        Parameters
        ----------
        x: Index to get unit network from data.

        Returns
        -------
        estimates: the beta and alpha
        """

        mu, var = self.__un[x].mean(), self.__un[x].var()
        al = ((1 - mu) / var - 1 / mu) * mu ** 2
        be = al * (1 / mu - 1)
        return al, be

    def __unit_network(self):
        """
        Construct unit wise regulatory networks.

        """

        exp1 = get_scaled_values(self.x)
        w1 = self.w.select_dtypes(include=["number", "float", "int"]).copy()
        exp_set = exp1.select_dtypes(include=["number", "float", "int"]).copy().T

        self.__un = exp1[["name"]]
        self.__ud = exp1[["name"]]
        for i in range(w1.shape[1]):
            x = np.array(w1.iloc[:, i]).reshape((-1, 1))

            lm_g = []
            lm_d = []
            for t in range(exp_set.shape[1]):
                y = np.array(exp_set.iloc[:, t]).reshape((-1, 1))
                lm_r = LinearRegression(n_jobs=-1, fit_intercept=True, copy_X=False)
                lm_r.fit(x, y)
                lm_g.append(lm_r.score(x, y))
                lm_d.append(float(lm_r.coef_))

            lm_g = pd.DataFrame(lm_g, columns=[w1.columns.values[i]])
            lm_d = pd.DataFrame(lm_d, columns=[w1.columns.values[i]])

            self.__un = pd.concat([self.__un, lm_g], axis=1)
            self.__ud = pd.concat([self.__ud, lm_d], axis=1)
        return

    def __score_class(self):
        """
        Call interactions {"P", "B} and
        counts them for global prediction.

        """

        nn = get_col_names(self.__un)
        for i in nn:
            a, b = self.__beta_params(i)
            self.__un[i] = 1 - beta.cdf(self.__un[i], a, b)
            self.__un[i] = np.where(self.__un[i] <= self.p, "P", "B")

        ae_col = [col for col in self.__un if col.startswith("AE")]
        filter_ae = self.__un[ae_col]
        self.__un["Predicted"] = filter_ae.iloc[:].eq("P").sum(axis=1).to_frame()
        return

    def fit(self):
        self.__unit_network()
        self.__score_class()
        return

    def save_net(self, s, i, nn):
        """
        Save the unit network with the directions.

        Parameters
        ----------
        i: Gene set name
        s: Boot number
        nn: Input name

        """
        self.__un.to_csv("N{}_{}_{}".format(s, i, nn), index=False)
        self.__ud.to_csv("D{}_{}_{}".format(s, i, nn), index=False)
        return

    def get_predicted(self, t):

        if t == 1:
            op = self.__un
        elif t == 0:
            op = self.__un[[col for col in self.__un.columns if col in ["name", "Predicted"]]]

        else:
            op = 0
        return op


