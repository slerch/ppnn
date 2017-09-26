"""
EMOS Network implementation in theano
"""

#Imports
import theano
import theano.tensor as T
import numpy as np


class EMOS_Network(object):
    def __init__(self):
        """
        This function is called once an object of this class is created.
        """
        # Before we start with the network, let's define
        # the learning rate as an input so we can vary it
        lr = T.fscalar('lr')
        
        # First let's define the input to the network
        # This is the ensemble mean (meanx), 
        # the ensemble stadnard deviation (stdx) and
        # the corresponding observation (target)
        # In theano we use tensors to describe these variables.
        # T.fvector allocates a float32 1D vector
        meanx = T.fvector('meanx')   # The name helps with debugging
        stdx = T.fvector('stdx')
        target = T.fvector('target')
        
        # Next we allocate the weights (a, b, c, d) as shared
        # variables and initialize some value for them.
        # For now we will just draw a random variable from N(0, 1)
        a = theano.shared(np.random.randn(), 'a')
        b = theano.shared(np.random.randn(), 'b')
        c = theano.shared(np.random.randn(), 'c')
        d = theano.shared(np.random.randn(), 'd')
        
        # Now that we have the input and the weights, 
        # we can set up the network.
        mu = a + meanx * b
        sigma = c + stdx * d
        out = T.stack([mu, sigma], axis=1)

        mean_crps = crps_cost_function(target, out, theano=True)

        # Now compute the gradients of the cost function 
        # with respect to the four weights/parameters
        params = [a, b, c, d]   # Let's put them in a list for convenience
        gradients = theano.tensor.grad(
            crps_cost_function(target, out, theano=True), 
            params,
            )
        
        # For gradient descent we now need to subtract the gradients
        # from our parameters to minimize the cost function
        # In theano we want to define a list of tuples containing
        # the old parameter and the updated parameter.
        updates = [(p, p - lr * g) for p, g in zip(params, gradients)]
        
        # So far no actual computations have been done. Now we will
        # define a Theano function, which takes input, does some 
        # calculations and returns some output. In our case, we use 
        # meanx, stdx and the target as an input plus the required 
        # learning rate and return the mean_crps
        # as an output. Then we tell the function to apply the update
        # every time it is called. This is the training
        self.train = theano.function([meanx, stdx, target, lr], 
                                     mean_crps, updates=updates)
        # Furthermore, we define a method for simply making a prediction
        # and returning the predicted values of mu and sigma
        # along with the mean_crps without updating the parameters
        self.predict = theano.function([meanx, stdx, target],
                                       [mu, sigma, mean_crps, a, b, c, d])

    def fit(self, meanx, stdx, target, epochs_max, validation_data, lr=0.1, 
            early_stopping_delta=None, verbose=0):
        """Function for fitting. 

        Similar to Keras fit method.
        Should also take batch_size, epochs. Shuffle
        """

        # Set up the early stopping
        # For now this is implemented to stop when the average train loss
        # stops decreasing by delta.
        # Average is over the last 5 epochs.
        tmp_cost_list = [1e99] * 5   # Start with a really high cost
        mean_cost_old = np.mean(tmp_cost_list[-5:])

        for i in range(epochs_max):
            # Train the model with gradient descent
            cost = self.train(meanx, stdx, target, lr)

            # Check if early stopping
            if early_stopping_delta is not None:
                tmp_cost_list.append(cost)
                mean_cost = np.mean(tmp_cost_list[-5:])
                if mean_cost_old - mean_cost < early_stopping_delta:
                    if verbose == 1: print('Stop training at step %i' % i)
                    break
                mean_cost_old = mean_cost

            if i%50 == 0 and verbose ==1: 
                print('Step %i; mean_crps = %.3f' % (i, cost))

        # Get the validation score
        crps_valid = self.predict(validation_data[0], validation_data[1], 
                                  validation_data[2])[2]

        return cost, crps_valid



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
    var = T.sqr(sigma)
    # The following three variables are just for convenience
    loc = (target - mu) / T.sqrt(var)
    phi = 1.0 / np.sqrt(2.0 * np.pi) * T.exp(-T.square(loc) / 2.0)
    Phi = 0.5 * (1.0 + T.erf(loc / np.sqrt(2.0)))
    # First we will compute the crps for each input/target pair
    crps =  T.sqrt(var) * (loc * (2. * Phi - 1.) + 2 * phi - 1. / np.sqrt(np.pi))
    # Then we take the mean. The cost is now a scalar
    return T.mean(crps)

def crps_cost_function_var(target, pred, theano=False):
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
    var = pred[:, 1]
    if not theano:   # Ugly workaround for different tensor allocation in keras and theano
        target = target[:, 0]   # Need to also get rid of axis 1 to match!

    # To stop sigma from becoming negative we first have to 
    # convert it the the variance and then take the square
    # root again. 
 
    # The following three variables are just for convenience
    loc = (target - mu) / T.sqrt(var)
    phi = 1.0 / np.sqrt(2.0 * np.pi) * T.exp(-T.square(loc) / 2.0)
    Phi = 0.5 * (1.0 + T.erf(loc / np.sqrt(2.0)))
    # First we will compute the crps for each input/target pair
    crps =  T.sqrt(var) * (loc * (2. * Phi - 1.) + 2 * phi - 1. / np.sqrt(np.pi))
    # Then we take the mean. The cost is now a scalar
    return T.mean(crps)