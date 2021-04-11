"""
Name: ACSNI-database
Author: Chinedu A. Anene, Phd
"""

import pandas as pd
from sklearn.linear_model import LinearRegression,LogisticRegression
import numpy as np
import pickle
import sys


class AcsniResults:
  """
  Class for ACSNI results
  """

  def __init__(self, co, w, n, p, d, run_info):
      self.co = co
      self.w = w
      self.n = n
      self.p = p
      self.d = d
      self.__lm_w = None
      self.__nm = None
      self.run_info = run_info

  def __str__(self):
      n_boot = self.p["Boot_Count"][1]
      n_modules = int((len(self.w.columns)-1)/n_boot)
      return "ACSNI result with {} modules over {} bootstraps.".format(n_modules, n_boot)

  def __repr__(self):
      return "\n" + self.__str__()

  def get_co(self):
      return self.co

  def get_w(self):
      return self.w

  def get_n(self):
      return self.n

  def get_p(self):
      return self.p

  def get_d(self):
      return self.d

  def get_run_info(self):
      return self.run_info

  def get_global_prediction(self):
      gr = self.p.copy()
      gr["name"] = gr.index
      gr = pd.concat([self.p["name"],
                      self.p.select_dtypes(include=["number", "float", "int"]).copy()],
                     axis=1)
      return gr

  def get_n_boot_units(self):
      boot = self.p["Boot_Count"][1]
      mod = int((len(self.w.columns)-1)/boot)
      return boot, mod

  def get_associated_module(self, target):
      names_tar = target.loc[:, target.columns != "ID"].columns[0]
      tar = self.get_w()
      tar = pd.merge(tar, target, left_on="ID", right_on="ID", how="left")
      return tar, names_tar

  def get_sub_process_network(self, c):
      n_p = self.p.loc[self.p["Sum_stat"] >= c]
      n_p = n_p.select_dtypes(exclude=["number", "float", "int"]).copy()
      n_p = n_p[[col for col in n_p if col.startswith("AE")]]
      n_p["name"] = n_p.index
      n_p = n_p.melt(id_vars="name", var_name="sub", value_name="Predicted")
      n_p = n_p.loc[n_p["Predicted"] == "P"]
      n_p["Set"] = "w" + n_p["sub"].astype(str)
      n_p = n_p[["name", "sub", "Set"]]
      return n_p

  def add_sub_process_direction(self, c):
      n_p = self.get_sub_process_network(c)
      n_p["match"] = n_p["name"] + n_p["Set"]
      self.d["name"] = self.d.index
      d_p = self.d.melt(id_vars="name", var_name="sub", value_name="Direction")
      d_p["Set"] = "w" + d_p["sub"].astype(str)
      d_p["match"] = d_p["name"] + d_p["Set"]
      d_p = d_p.loc[d_p["match"].isin(n_p["match"])]
      d_p = d_p[["Direction", "match"]]
      n_p = pd.merge(n_p, d_p, how="left", on="match")
      return n_p[["name", "sub", "Direction"]]

  def get_explained_variation(self, t, c):
      """
      Get level of interaction between phenotype and subprocesses.

      """

      self.__nm = t.loc[:, t.columns != "ID"].columns[0]

      weights = self.w.select_dtypes(include=["number", "float", "int"]).columns
      weights = [cc for cc in weights if cc.startswith("AE")]
      lm_w = pd.DataFrame({"sub": weights})

      self.w = pd.merge(self.w, t, left_on="ID", right_on="ID", how="left")

      lm_g = []
      for i in weights:
          x = np.array(self.w[i]).reshape(-1, 1)
          y = self.w[[self.__nm]].values.ravel()

          if c == "character":
              lm_r = LogisticRegression(n_jobs=-1)
          elif c == "numeric":
              lm_r = LinearRegression(n_jobs=-1)
          else:
              sys.exit("Invalid value set for parameter -c. "
                       "Use numeric or character for the type of phenotype")

          lm_r.fit(x, y)
          lm_r.score(x, y)
          lm_g.append(lm_r.score(x, y))
      lm_g = pd.DataFrame(lm_g)
      self.__lm_w = pd.concat([lm_w, lm_g], axis=1)
      return

  def get_summary_stat_phenotype(self, t, c):
      """
      Phenotype association statistics
      Parameters
      ----------
      t: phenotype file
      c: Type of variable

      Returns
      -------
      results and name

      """

      self.get_explained_variation(t, c)
      qq = self.__lm_w[[0]].quantile([.25, .75])
      q25 = qq.loc[0.25, 0]
      q75 = qq.loc[0.75, 0]

      m = float(self.__lm_w[[0]].mean())
      sd = float(self.__lm_w[[0]].std())
      print("Statistics of variations in subprocesses explained by {}".format(self.__nm))
      print("q25 {} \n q75 {} \n mean {} \n std {}".format(q25, q75, m, sd))
      self.__lm_w["Association"] = np.where(self.__lm_w[0] >= q75, "Strong", "Weak")
      return self.__lm_w, self.__nm


class ACSNIDeriveResults:
  """
  Class for ACSNI Derive results
  """

  def __init__(self, ac, n, d, fd, co, run_info):
      self.n = n
      self.d = d
      self.ac = ac
      self.fd = fd
      self.co = co
      self.run_info = run_info

  def __str__(self):
      return "ACSNI Derive raw results"

  def __repr__(self):
      return "\n" + self.__str__()

  def get_run_info(self):
      return self.run_info

  def get_ac(self):
      return self.ac

  def get_n(self):
      return self.n

  def get_d(self):
      return self.d

  def get_fd(self):
      return self.fd

  def get_co(self):
      return self.co

def save_acsni(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

#dbst = load_acsni("dbsAA_.csv.ptl")
