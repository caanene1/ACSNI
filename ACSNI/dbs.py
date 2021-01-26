"""
Name: ACSNI-database
Author: Chinedu A. Anene, Phd
"""

import pandas as pd
import pickle


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
      gr = pd.concat([self.p["gene"],
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
      n_p = n_p.melt(id_vars="gene", var_name="sub", value_name="Predicted")
      n_p = n_p.loc[n_p["Predicted"] == "P"]
      n_p["Set"] = "w" + n_p["sub"].astype(str)
      n_p = n_p[["gene", "sub", "Set"]]
      return n_p

  def add_sub_process_direction(self, c):
      n_p = self.get_sub_process_network(c)
      n_p["match"] = n_p["gene"] + n_p["Set"]
      d_p = self.d.melt(id_vars="gene", var_name="sub", value_name="Direction")
      d_p["Set"] = "w" + d_p["sub"].astype(str)
      d_p["match"] = d_p["gene"] + d_p["Set"]
      d_p = d_p.loc[d_p["match"].isin(n_p["match"])]
      d_p = d_p[["Direction", "match"]]
      n_p = pd.merge(n_p, d_p, how="left", on="match")
      return n_p


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
