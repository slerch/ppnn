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


def approx_crps_cat(y_true, y_pred):
    # [sample, cat]
    cum_obs = K.cumsum(y_true, axis=1)
    cum_preds = K.cumsum(y_pred, axis=1)
    approx_crps = K.sum(K.square(cum_obs - cum_preds), axis=1) * bin_width
    approx_crps = K.mean(approx_crps)
    return approx_crps


def maybe_correct_cat_crps(preds, targets, bin_edges):
    """CRPS for categorical predictions. Not sure if correct

    """
    # pdb.set_trace()
    # Convert input arrays
    preds = np.array(np.atleast_2d(preds), dtype='float')
    targets = np.array(np.atleast_1d(targets), dtype='float')

    # preds [sample, bins]
    # Find insert index
    mat_bins = np.repeat(np.atleast_2d(bin_edges), targets.shape[0], axis=0)
    b = mat_bins.T - targets
    b[b < 0] = 999
    insert_idxs = np.argmin(b, axis=0)

    # Insert
    ins_bin_edges = np.array([np.insert(np.array(bin_edges, dtype=float),
                                        insert_idxs[i], targets[i])
                              for i in range(targets.shape[0])])
    ins_preds = np.array(
        [np.insert(preds[i], insert_idxs[i], preds[i, insert_idxs[i] - 1])
         for i in range(targets.shape[0])])

    # Get obs
    bin_obs = np.array([(ins_bin_edges[i, :-1] <= targets[i]) &
                        (ins_bin_edges[i, 1:] > targets[i])
                        for i in range(targets.shape[0])], dtype=int)

    # Cumsum with weights
    ins_preds *= np.diff(ins_bin_edges, axis=1)
    cum_bin_obs = np.cumsum(bin_obs, axis=1)
    cum_probs = np.cumsum(ins_preds, axis=1)
    cum_probs = (cum_probs.T / cum_probs[:, -1]).T

    # Get adjusted preds
    adj_cum_probs = np.concatenate(
        (np.zeros((cum_probs.shape[0], 1)), cum_probs), axis=1)
    adj_cum_probs = (adj_cum_probs[:, :-1] + adj_cum_probs[:, 1:]) / 2

    # Compute CRPS
    crps = np.mean(np.sum(((adj_cum_probs - cum_bin_obs) ** 2) *
                          np.diff(ins_bin_edges, axis=1), axis=1))
    return crps
