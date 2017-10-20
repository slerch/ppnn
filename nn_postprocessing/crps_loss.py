"""
Definition of CRPS loss function.
"""
import keras
import keras.backend as K
import theano.tensor as T
import numpy as np
if keras.backend.backend() == 'tensorflow':
    from tensorflow import erf
else:
    from theano.tensor import erf

def crps_cost_function(target, pred, theano=False):
    """Computes the mean CRPS loss function.
    
    Code inspired by Kai Polsterer (HITS)
    Inputs
    ------
        pred = [mu, sigma] : Theano tensors containing means and stds as np array
        target: Theano tensor containing target
    Outputs
    -------
        mean_crps : Scalar with mean CRPS over all samples
    """
    # Split input
    mu = pred[:, 0]
    sigma = pred[:, 1]
    if not theano:   # Ugly workaround for different tensor allocation in keras and theano
        target = target[:, 0]   # Need to also get rid of axis 1 to match!

    # To stop sigma from becoming negative we first have to 
    # convert it the the variance and then take the square
    # root again. 
    var = K.square(sigma)
    # The following three variables are just for convenience
    loc = (target - mu) / K.sqrt(var)
    phi = 1.0 / np.sqrt(2.0 * np.pi) * K.exp(-K.square(loc) / 2.0)
    Phi = 0.5 * (1.0 + erf(loc / np.sqrt(2.0)))
    # First we will compute the crps for each input/target pair
    crps =  K.sqrt(var) * (loc * (2. * Phi - 1.) + 2 * phi - 1. / np.sqrt(np.pi))
    # Then we take the mean. The cost is now a scalar
    return K.mean(crps)


def crps_cost_function_seq(target, pred):
    """TODO
    Inputs
    ------
        pred = [mu, sigma] : Theano tensors containing means and stds as np array
        target: Theano tensor containing target
    Outputs
    -------
        mean_crps : Scalar with mean CRPS over all samples
    """
    # Split input
    mu = pred[:, :, 0]
    sigma = pred[:, :, 1]
    
    tar = target[:, :, 0]
    # [sample, time_step]

    # To stop sigma from becoming negative we first have to 
    # convert it the the variance and then take the square
    # root again. 
    var = K.square(sigma)
    # The following three variables are just for convenience
    loc = (tar - mu) / K.sqrt(var)
    phi = 1.0 / np.sqrt(2.0 * np.pi) * K.exp(-K.square(loc) / 2.0)
    Phi = 0.5 * (1.0 + erf(loc / np.sqrt(2.0)))
    # First we will compute the crps for each input/target pair
    crps =  K.sqrt(var) * (loc * (2. * Phi - 1.) + 2 * phi - 1. / np.sqrt(np.pi))

    # Here we do not take the mean because we want keras to be able to apply
    # weights afterwards!
    return crps