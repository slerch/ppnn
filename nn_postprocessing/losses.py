"""
Definition of CRPS loss function.
"""
import keras
import keras.backend as K
import numpy as np
# Import erf depending on whether we use the theano or tensorflow backend
if keras.backend.backend() == 'tensorflow':
    from tensorflow import erf
else:
    from theano.tensor import erf


def crps_cost_function(y_true, y_pred, theano=False):
    """Compute the CRPS cost function for a normal distribution defined by
    the mean and standard deviation.

    Code inspired by Kai Polsterer (HITS).

    Args:
        y_true: True values
        y_pred: Tensor containing predictions: [mean, std]
        theano: Set to true if using this with pure theano.

    Returns:
        mean_crps: Scalar with mean CRPS over batch
    """

    # Split input
    mu = y_pred[:, 0]
    sigma = y_pred[:, 1]
    # Ugly workaround for different tensor allocation in keras and theano
    if not theano:
        y_true = y_true[:, 0]   # Need to also get rid of axis 1 to match!

    # To stop sigma from becoming negative we first have to 
    # convert it the the variance and then take the square
    # root again. 
    var = K.square(sigma)
    # The following three variables are just for convenience
    loc = (y_true - mu) / K.sqrt(var)
    phi = 1.0 / np.sqrt(2.0 * np.pi) * K.exp(-K.square(loc) / 2.0)
    Phi = 0.5 * (1.0 + erf(loc / np.sqrt(2.0)))
    # First we will compute the crps for each input/target pair
    crps =  K.sqrt(var) * (loc * (2. * Phi - 1.) + 2 * phi - 1. / np.sqrt(np.pi))
    # Then we take the mean. The cost is now a scalar
    return K.mean(crps)


def crps_cost_function_seq(y_true, y_pred):
    """Version of CRPS const function for sequence predictions.
    Here the input tensors have dimensions [sample, time_step].
    The output has the same dimensions so that keras can apply weights
    afterwards for missing data.

    Args:
        y_true: True values with dimensions [sample, time_step, 1]
        y_pred: Predictions with dimensions [sample, time_step, [mean, std]]

    Returns:
        crps: CRPS with dimensions [sample, time_step]
    """
    # Split input
    mu = y_pred[:, :, 0]
    sigma = y_pred[:, :, 1]
    
    tar = y_true[:, :, 0]
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


def approx_crps_cat(bin_width):
    """Wrapper to pass bin_width as an argument to the loss function.

    Args:
        bin_width: width of categorical bins

    Returns:
        loss_function: approximate crps_loss function with bin_width specified
    """
    def loss(y_true, y_pred):
        """Approximate CRPS function for categorical output.

        Args:
            y_true: One-hot-encoded output
            y_pred: Probability for each bin

        Returns:
            approx_crps: Approximate mean CRPS value for batch
        """
        # [sample, cat]
        cum_obs = K.cumsum(y_true, axis=1)
        cum_preds = K.cumsum(y_pred, axis=1)
        approx_crps = K.sum(K.square(cum_obs - cum_preds), axis=1) * bin_width
        return K.mean(approx_crps)
    return loss
