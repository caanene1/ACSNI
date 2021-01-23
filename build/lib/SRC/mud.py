"""
Name: ACSNI-model
Author: Chinedu A. Anene, Phd
"""

from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras import regularizers, optimizers
import numpy as np
from pandas import DataFrame as Dff
from SRC.dat import name_generator


def aer(nc, nr, lp):
    """
    Build model structure

    Parameters
    ----------
    nc: Number of columns
    nr: Number of rows
    lp: The proportion of the gene-set to use as the layer size

    Returns
    -------
    auto_encoder: Data Specific model
    """
    print("Compiling model structure")
    a = nr * lp // 100
    b = a // 2

    input_layer = Input(shape=(nc,), name="input")
    encoder = Dense(a, activation="relu", kernel_initializer="glorot_uniform",
                    activity_regularizer=regularizers.l1_l2(1e-16, 1e-9),
                    name="enco1")(input_layer)

    encoder = Dense(b, activation="relu", name="code")(encoder)

    decoder = Dense(a, activation="sigmoid", name="deco1")(encoder)
    decoder = Dense(nc, activation="sigmoid", name="output")(decoder)

    auto_encoder = Model(inputs=input_layer, outputs=decoder)
    auto_encoder.compile(optimizer= optimizers.RMSprop(learning_rate=1e-3),
                         loss='mean_squared_error', metrics=['mse'])
    # auto_encoder.summary()
    return auto_encoder

def mud_output(x_train, w_nn, lp):
    """
    Model training and output.

    Parameters
    ----------
    x_train: Training data set
    w_nn: Optional weights for the genes
    lp: scaling percentage

    Returns
    -------
    final_model: The best model
    weights: Estimated weights for the samples
    code: Pathway representation
    """
    np.random.seed(100)
    run_id = name_generator(6)
    nam = run_id  + "_model.h5"

    final_model = aer(nc=x_train.shape[1], nr=x_train.shape[0], lp=lp)

    check = ModelCheckpoint(filepath=nam, verbose=0, save_best_only=True)

    final_model.fit(x_train, x_train, sample_weight=w_nn,
                    epochs=3000, shuffle=True, batch_size=15,
                    validation_data=(x_train, x_train), callbacks=[check],
                    verbose=0)

    final_model = load_model(nam)
    weights = final_model.get_layer("enco1").get_weights()[0]

    code = final_model.get_layer("deco1").get_weights()[0]
    code = Dff(code)
    code["run"] = nam
    code.to_csv("code_{}.csv".format(run_id), index=False)
    return weights
