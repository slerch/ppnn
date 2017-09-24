"""
Utility functions for PPNN.
The functions are still a little confusingly named and 
structured at the moment.

Author: Stephan Rasp
"""

from scipy.stats import norm
import numpy as np
from netCDF4 import num2date, Dataset
from emos_network_theano import EMOS_Network
import timeit


def load_nc_data(fn, utc=0):
    """
    Returns the full dataset from the netCDF file for a given valid time.
    So: utc = 0 or 12
    tobs, tfc, dates
    """
    rg = Dataset(fn)

    # Load the data as (masked) numpy arrays
    tobs = rg.variables['t2m_obs'][:]
    tfc = rg.variables['t2m_fc'][:]
    dates = num2date(rg.variables['time'][:],
                     units='seconds since 1970-01-01 00:00 UTC')

    # For now, convert data to C
    idx = np.where(np.mean(tfc, axis=(1, 2)) > 100)[0][0]
    tfc[idx:] = tfc[idx:] - 273.15

    # Compute hours
    hours = np.array([d.hour for d in list(dates)])

    # Get data for given valid time
    tobs = tobs[hours == 0]
    tfc = tfc[hours == 0]
    dates = dates[hours == 0]

    return tobs, tfc, dates


def get_train_test_data(tobs_full, tfc_full, date_idx, window_size=25, fclt=48):
    """
    Returnes the prepared and normalized training and test data.
    Training data: tobs and tfc for the rolling window
    Test data: tobs and tfc for the to be predicted date.
    """
    # Get the data from the full data set
    tobs_train, tfc_train = get_rolling_slice(tobs_full, tfc_full, date_idx, 
                                            window_size, fclt)
    tobs_test, tfc_test = (tobs_full[date_idx], tfc_full[date_idx])

    # Compress the data and remove nans
    tobs_train, tfc_mean_train, tfc_std_train = prep_data(tobs_train, tfc_train)
    tobs_test, tfc_mean_test, tfc_std_test = prep_data(tobs_test, tfc_test)

    # Scale the input features
    tfc_mean_mean = tfc_mean_train.mean()
    tfc_mean_std = tfc_mean_train.std()
    tfc_std_mean = tfc_std_train.mean()
    tfc_std_std = tfc_std_train.std()

    tfc_mean_train = (tfc_mean_train - tfc_mean_mean) / tfc_mean_std
    tfc_mean_test = (tfc_mean_test - tfc_mean_mean) / tfc_mean_std
    tfc_std_train = (tfc_std_train - tfc_std_mean) / tfc_std_std
    tfc_std_test = (tfc_std_test - tfc_std_mean) / tfc_std_std

    return (tfc_mean_train, tfc_std_train, tobs_train, 
            tfc_mean_test, tfc_std_test, tobs_test)



def get_rolling_slice(tobs_full, tfc_full, date_idx, window_size=25, fclt=48):
    """
    Return the forecast and observation data from the 
    previous *window_size* days. So if date_idx=10 and 
    window_size=3, it would get the data for indices 7, 8, 9.
    Nope, also have to go back the forecast lead time.
    """
    fclt_didx = int(fclt / 24)
    
    # Get the correct indices
    idx_start = date_idx - window_size - fclt_didx
    idx_stop = date_idx - fclt_didx
    
    # Get the slice for the indices
    tobs_roll = tobs_full[idx_start:idx_stop]
    tfc_roll = tfc_full[idx_start:idx_stop]
    
    return tobs_roll, tfc_roll


def get_data_slice(rg, month, utc=0):
    """
    Get data for one month.

    Not currently used!
    """
    # Get array of datetime objects
    dates = num2date(rg.variables['time'][:], 
                     units='seconds since 1970-01-01 00:00 UTC')
    # Extract months and hours
    months = np.array([d.month for d in list(dates)])
    hours = np.array([d.hour for d in list(dates)])
    
    # for now I need to include the Kelvin fix
    tfc = rg.variables['t2m_fc'][:]
    idx = np.where(np.mean(tfc, axis=(1, 2)) > 100)[0][0]
    tfc[idx:] = tfc[idx:] - 273.15
    
    # Extract the requested data
    tobs = rg.variables['t2m_obs'][(months == month) & (hours == utc)]
    tfc = tfc[(months == month) & (hours == utc)]
    
    return tobs, tfc


def prep_data(tobs, tfc, verbose=False):
    """
    Prepare the data as input for Network.
    """
    ax = 0 if tobs.ndim == 1 else 1
    # Compute mean and std and convert to float32
    tfc_mean = np.mean(np.asarray(tfc, dtype='float32'), axis=ax)
    tfc_std = np.std(np.asarray(tfc, dtype='float32'), axis=ax, ddof=1)
    tobs = np.asarray(tobs, dtype='float32')
    
    # Flatten
    tobs = np.ravel(tobs)
    tfc_mean = np.ravel(tfc_mean)
    tfc_std = np.ravel(tfc_std)
    
    # Remove NaNs
    mask = np.isfinite(tobs)
    if verbose:
        print('NaNs / Full = %i / %i' % (np.sum(~mask), tobs.shape[0]))
    tobs = tobs[mask]
    tfc_mean = tfc_mean[mask]
    tfc_std = tfc_std[mask]
    
    return tobs, tfc_mean, tfc_std


def crps_normal(mu, sigma, y):
    """
    Compute CRPS for a Gaussian distribution. 
    """
    loc = (y - mu) / sigma
    crps = sigma * (loc * (2 * norm.cdf(loc) - 1) + 
                    2 * norm.pdf(loc) - 1. / np.sqrt(np.pi))
    return crps


def loop_over_days(tobs_full, tfc_full, date_idx_start, date_idx_stop, 
                   window_size, fclt, epochs_max, early_stopping_delta=None,
                   lr=0.1, reinit_model=False, verbose=0):
    """Function to loop over days with Theano EMOS_Network model.

    Parameters
    ----------
    tobs_full : numpy array [time, station]
        Output of load_nc_data
    tfc_full : numpy array [time, member, station]
        Output of load_nc_data
    date_idx_start/stop : int
        Start and stop index for loop
    window_size : int
        How many days for rolling training period
    fclt : int
        Forecast lead time in hours
    epochs_max : int
        How many times to fit to the entire training set
    early_stopping_delta : float
        Minimum improvement in mean train CRPS to keep going
    lr : float
        Learning rate
    reinit_model : bool
        If True, model weights are reinitialized for each day.
        If False, model weights are kept and just updated
    verbose : int
        0 or 1. If 1, print additional output
    Returns
    -------
    train_crps_list : list
        List with training CRPS for each day
    valid_crps_list : list
        List with validation/prediction CRPS for each day
    """

    # Make sure learning rate is 32 bit float
    lr = np.asarray(lr, dtype='float32')

    # Set up model initially
    model = EMOS_Network()

    # Allocate lists to write results
    train_crps_list = []
    valid_crps_list = []

    # Start timer and loop 
    time_start = timeit.default_timer()
    for date_idx in range(date_idx_start, date_idx_stop + 1):
        if date_idx % 20 == 0:
            print(date_idx)

        # Get data slice
        tfc_mean_train, tfc_std_train, tobs_train, \
            tfc_mean_test, tfc_std_test, tobs_test = \
            get_train_test_data(tobs_full, tfc_full, date_idx, 
                                window_size, fclt)

        # Reinitialize model if requested
        if reinit_model:
            model = EMOS_Network()

        # Fit model
        train_crps, valid_crps = model.fit(
            tfc_mean_train, tfc_std_train, tobs_train, 
            epochs_max, 
            (tfc_mean_test, tfc_std_test, tobs_test), 
            lr=lr, 
            early_stopping_delta=early_stopping_delta
            )

        # Write output
        train_crps_list.append(train_crps)
        valid_crps_list.append(valid_crps)

    # Stop timer 
    time_stop = timeit.default_timer()
    print('Time: %.2f s' % (time_stop - time_start))

    return train_crps_list, valid_crps_list




