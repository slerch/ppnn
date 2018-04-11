import numpy as np
import sys, os, pdb, pickle, time
from .utils import *
from .losses import *
from .keras_models import * 
from .aux_dict import *
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib import animation
import seaborn as sns
from tqdm import tqdm_notebook as tqdm
from collections import OrderedDict
from IPython.display import HTML
from keras.utils.generic_utils import get_custom_objects
metrics_dict = dict([(f.__name__, f) for f in [crps_cost_function]])
get_custom_objects().update(metrics_dict)
from timeit import default_timer
from keras.callbacks import EarlyStopping
import tensorflow as tf

def limit_mem():
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    keras.backend.tensorflow_backend.set_session(tf.Session(config=config))
limit_mem()