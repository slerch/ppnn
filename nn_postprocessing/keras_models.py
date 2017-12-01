"""Definition of keras models used in the Jupyter notebooks.

"""

from losses import crps_cost_function

import keras
from keras.layers import Input, Dense, merge, Embedding, Flatten, Concatenate
from keras.models import Model
import keras.backend as K
from keras.callbacks import EarlyStopping


def build_EMOS_network_keras(compile=False, optimizer='sgd', lr=0.1):
    """Build (and maybe compile) EMOS network in keras.

    Args:
        compile: If true, compile model
        optimizer: String of keras optimizer
        lr: learning rate

    Returns:
        model: Keras model
    """
    mean_in = Input(shape=(1,))
    std_in = Input(shape=(1,))
    mean_out = Dense(1, activation='linear')(mean_in)
    std_out = Dense(1, activation='linear')(std_in)
    x = keras.layers.concatenate([mean_out, std_out], axis=1)
    model = Model(inputs=[mean_in, std_in], outputs=x)

    if compile:
        opt = keras.optimizers.__dict__[optimizer](lr=lr)
        model.compile(optimizer=opt, loss=crps_cost_function)
    return model


def build_fc_model(n_features, n_outputs, compile=False, optimizer='adam',
                   lr=0.1, loss=crps_cost_function):
    """Build (and compile) fully connected linear network.

    Args:
        n_features: Number of features
        n_outputs: Number of outputs
        compile: If true, compile model
        optimizer: Name of optimizer
        lr: learning rate
        loss: loss function

    Returns:
        model: Keras model
    """

    inp = Input(shape=(n_features,))
    x = Dense(n_outputs, activation='linear')(inp)
    model = Model(inputs=inp, outputs=x)

    if compile:
        opt = keras.optimizers.__dict__[optimizer](lr=lr)
        model.compile(optimizer=opt, loss=loss)
    return model


def build_hidden_model(n_features, n_outputs, hidden_nodes, compile=False,
                       optimizer='adam', lr=0.01, loss=crps_cost_function,
                       activation='relu'):
    """Build (and compile) a neural net with hidden layers

    Args:
        n_features: Number of features
        n_outputs: Number of outputs
        hidden_nodes: int or list of hidden nodes
        compile: If true, compile model
        optimizer: Name of optimizer
        lr: learning rate
        loss: loss function

    Returns:
        model: Keras model
    """
    if type(hidden_nodes) is not list:
        hidden_nodes = [hidden_nodes]
    inp = Input(shape=(n_features,))
    x = Dense(hidden_nodes[0], activation=activation)(inp)
    if len(hidden_nodes) > 1:
        for h in hidden_nodes[1:]:
            x = Dense(h, activation=activation)(x)
    x = Dense(n_outputs, activation='linear')(x)
    model = Model(inputs=inp, outputs=x)

    if compile:
        opt = keras.optimizers.__dict__[optimizer](lr=lr)
        model.compile(optimizer=opt, loss=loss)
    return model


def build_emb_model(n_features, n_outputs, hidden_nodes, emb_size, max_id,
                    compile=False, optimizer='adam', lr=0.01,
                    loss=crps_cost_function):
    """

    Args:
        n_features: Number of features
        n_outputs: Number of outputs
        hidden_nodes: int or list of hidden nodes
        emb_size: Embedding size
        max_id: Max embedding ID
        compile: If true, compile model
        optimizer: Name of optimizer
        lr: learning rate
        loss: loss function

    Returns:
        model: Keras model
    """
    if type(hidden_nodes) is not list:
        hidden_nodes = [hidden_nodes]

    features_in = Input(shape=(n_features,))
    id_in = Input(shape=(1,))
    emb = Embedding(max_id + 1, emb_size)(id_in)
    emb = Flatten()(emb)
    x = Concatenate()([features_in, emb])
    for h in hidden_nodes:
        x = Dense(h, activation='relu')(x)
    x = Dense(n_outputs, activation='linear')(x)
    model = Model(inputs=[features_in, id_in], outputs=x)

    if compile:
        opt = keras.optimizers.__dict__[optimizer](lr=lr)
        model.compile(optimizer=opt, loss=loss)
    return model