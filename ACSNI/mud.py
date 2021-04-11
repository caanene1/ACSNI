"""
Name: ACSNI-model
Author: Chinedu A. Anene, Phd
"""
from sklearn.decomposition import PCA, NMF
import warnings
from sklearn import metrics
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras import regularizers, optimizers
import numpy as np
from pandas import DataFrame as Dff
import pandas as pd
from ACSNI.dat import name_generator



class DimRed:
  """
  Class for quadruple dimension reduction.
  """

  def __init__(self, x, w, p):
      """
      Dimension reduction class.

      Parameters:
          x: Input matrix (np array)
          w: Weights to adjustment for ae
          p: latent dimension for ae
      """
      self.x = x
      self.w = w
      self.p = p
      self.pca = None
      self.nmf = None
      self.ae = None
      self.__reduced = None
      self.__pcanmf = None
      self.__median = None
      self.__ael1 = None
      self.__ael2 = None
      self.__a = None
      self.__r = 15
      self.__a1 = 0.04
      self.__a2 = 0.80
      self.__scorer = metrics.explained_variance_score
      self.__run_id = name_generator(6)

  def __str__(self):
      return "Quadruple dimension reduction class"

  def __repr__(self):
      return "\n" + self.__str__()

  def __get_score(self, model, y):
      """
      Determine level of explained variance

      """
      prediction = model.inverse_transform(model.transform(y))
      return self.__scorer(y, prediction)

  def l_med(self):
      self.__median = pd.DataFrame(np.median(self.x.T, axis=1))
      self.__median = self.__median.add_prefix('MEDIAN_' + self.__run_id + '_')
      return


  def lde(self):
      """
      Decompose with PCA and NMF

      """
      self.pca = PCA(n_components=0.95)
      self.pca.fit(self.x)
      pc_weights = pd.DataFrame(self.pca.components_.T)
      pc_weights = pc_weights.add_prefix('PCA_' + self.__run_id + '_')

      opti_rank = []
      #
      warnings.filterwarnings('ignore')
      for k in range(2, self.__r):
          nmf = NMF(n_components=k, max_iter=1000).fit(self.x)
          score_it = self.__get_score(nmf, self.x)
          opti_rank.append(score_it)
          if score_it >= 0.95:
              break


      self.nmf = NMF(n_components=len(opti_rank) + 1, max_iter=10000)
      self.nmf.fit(self.x)
      warnings.resetwarnings()
      #

      nmf_weights = pd.DataFrame(self.nmf.components_.T)
      nmf_weights = nmf_weights.add_prefix('NMF_' + self.__run_id + '_')

      self.__pcanmf = pd.concat([pc_weights, nmf_weights], axis=1)

      return

  def __de4ae(self, y):
      """
      Estimate optimal dimension for AE,
        based on Bahadur and Paffenroth 2020, IEEE

      """
      s_x = y.copy()

      for t in range(s_x.shape[0]):
          s_x.iloc[t, :] = np.sort(np.array(s_x.iloc[t, :]))[::-1]

      svp = np.sort(s_x.mean())[::-1]
      svp_sum = svp.sum()
      alg1 = sum(svp / svp_sum > self.__a1)
      alg2 = 0

      temp = (svp_sum * self.__a2) / 1
      temp2 = 0

      for i in range(len(svp)):
          temp2 += svp[i]
          alg2 += 1

          if temp2 >= temp:
              break
      return int((alg1 + alg2) / 2)

  def __aer(self, nc):
      """
      Build model structure

      """
      input_layer = Input(shape=(nc,), name="input")
      encoder = Dense(self.__a, activation="relu",
                      kernel_initializer="glorot_uniform",
                      activity_regularizer=regularizers.l1_l2(1e-16, 1e-9),
                      name="enco1")(input_layer)

      encoder = Dense(self.__a // 2, activation="relu", name="code")(encoder)

      decoder = Dense(self.__a, activation="sigmoid", name="deco1")(encoder)
      decoder = Dense(nc, activation="sigmoid", name="output")(decoder)

      self.ae = Model(inputs=input_layer, outputs=decoder)
      self.ae.compile(optimizer=optimizers.RMSprop(learning_rate=1e-3),
                           loss='mean_squared_error', metrics=['mse'])
      return

  def mud_ae(self):
      """
      Model training and output.

      """
      nam = self.__run_id + "_model.h5"

      if self.p == 0:
          self.__a = self.x.shape[0] * 50 // 100
          check1 = ModelCheckpoint(filepath="est_" + nam, verbose=0,
                                   save_best_only=True)

          self.__aer(nc=self.x.shape[1])
          self.ae.fit(self.x, self.x, sample_weight=self.w,
                       epochs=500, shuffle=True, batch_size=15,
                       validation_data=(self.x, self.x), callbacks=[check1],
                       verbose=0)

          estimate = load_model("est_" + nam)
          code_est = Dff(estimate.get_layer("enco1").get_weights()[0])
          self.__a = self.__de4ae(code_est)
          print("The optimal number of dimension is {}".format(self.__a))

      else:
          self.__a = self.p

      check = ModelCheckpoint(filepath=nam, verbose=0, save_best_only=True)

      self.__aer(nc=self.x.shape[1])
      self.ae.fit(self.x, self.x, sample_weight=self.w,
                epochs=3000, shuffle=True, batch_size=15,
                validation_data=(self.x, self.x), callbacks=[check],
                verbose=0)

      final = load_model(nam)
      self.__ael1 = Dff(final.get_layer("enco1").get_weights()[0])
      self.__ael1 = self.__ael1.add_prefix('AE_' + self.__run_id + '_')

      self.__ael2 = Dff(final.get_layer("code").get_weights()[0])
      self.__ael2["run"] = nam

      self.__ael2.to_csv("code_{}.csv".format(self.__run_id), index=False)
      return

  def fit(self):
      """
      Fit quadruple dimension reduction {Median, PCA, NMF, AE[DE]}

      """
      self.l_med()
      self.lde()
      self.mud_ae()
      self.__reduced = pd.concat([self.__ael1,
                                  self.__pcanmf,
                                  self.__median],
                                 axis=1)
      return

  def get_reduced(self):
      """
      Get reduced dimension

      """
      return self.__reduced

  def get_aede(self):
      return self.__a

  def add_reduced_row(self, y):
      self.__reduced["ID"] = y
      self.__pcanmf["ID"] = y
      self.__ael1["ID"] = y
      return

