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
from keras.callbacks import EarlyStopping
from datetime import datetime
import pdb


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

    # Compute hours
    hours = np.array([d.hour for d in list(dates)])

    # Get data for given valid time
    tobs = tobs[hours == utc]
    tfc = tfc[hours == utc]
    dates = dates[hours == utc]

    return tobs, tfc, dates


def get_train_test_data(tobs_full, tfc_full, date_idx, window_size=25, fclt=48,
                        subtract_std_mean=True, test_plus=1):
    """
    Returnes the prepared and normalized training and test data.
    Training data: tobs and tfc for the rolling window
    Test data: tobs and tfc for the to be predicted date.
    """
    # Get the data from the full data set
    tobs_train, tfc_train = get_rolling_slice(tobs_full, tfc_full, date_idx, 
                                            window_size, fclt)
    tobs_test, tfc_test = (tobs_full[date_idx:date_idx + test_plus], 
                           tfc_full[date_idx:date_idx + test_plus])

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
    if subtract_std_mean:
        tfc_std_train = (tfc_std_train - tfc_std_mean) / tfc_std_std
        tfc_std_test = (tfc_std_test - tfc_std_mean) / tfc_std_std
    else:
        tfc_std_train = tfc_std_train / tfc_std_std
        tfc_std_test = tfc_std_test / tfc_std_std

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

def return_date_idx(dates, year, month, day):
    return np.where(dates == datetime(year, month, day, 0, 0))[0][0]


def crps_normal(mu, sigma, y):
    """
    Compute CRPS for a Gaussian distribution. 
    """
    loc = (y - mu) / sigma
    crps = sigma * (loc * (2 * norm.cdf(loc) - 1) + 
                    2 * norm.pdf(loc) - 1. / np.sqrt(np.pi))
    return crps


def loop_over_days(model, tobs_full, tfc_full, date_idx_start, date_idx_stop, 
                   window_size, fclt, epochs_max, early_stopping_delta=None,
                   lr=0.1, reinit_model=False, verbose=0, model_type='keras'):
    """Function to loop over days with Theano EMOS_Network model.

    Parameters
    ----------
    model : model
        Model with fit method, either keras or theano
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
    model_type : str
        'keras' or 'theano'
    Returns
    -------
    train_crps_list : list
        List with training CRPS for each day
    valid_crps_list : list
        List with validation/prediction CRPS for each day
    """

    # Make sure learning rate is 32 bit float
    lr = np.asarray(lr, dtype='float32')

    # Allocate lists to write results
    train_crps_list = []
    valid_crps_list = []

    # Start timer and loop 
    time_start = timeit.default_timer()
    for date_idx in range(date_idx_start, date_idx_stop + 1):
        if date_idx % 100 == 0:
            print(date_idx)

        # Get data slice
        tfc_mean_train, tfc_std_train, tobs_train, \
            tfc_mean_test, tfc_std_test, tobs_test = \
            get_train_test_data(tobs_full, tfc_full, date_idx, 
                                window_size, fclt, subtract_std_mean=False)

        # Reinitialize model if requested
        # only works for theano model
        if reinit_model:
            assert model_type == 'thano', 'Reinitialize does not work with keras'
            model = EMOS_Network()

        # Fit model
        if model_type == 'EMOS_Network_theano':
            train_crps, valid_crps = model.fit(
                tfc_mean_train, tfc_std_train, tobs_train, 
                epochs_max, 
                (tfc_mean_test, tfc_std_test, tobs_test), 
                lr=lr, 
                early_stopping_delta=early_stopping_delta,
                verbose=verbose,
                )
        elif model_type == 'EMOS_Network_keras':
            es = EarlyStopping(monitor='loss', min_delta=early_stopping_delta, 
                               patience=2)
            batch_size=tfc_mean_train.shape[0]
            model.fit(
                [tfc_mean_train, tfc_std_train], tobs_train,
                epochs=epochs_max,
                batch_size=batch_size,
                verbose=verbose,
                callbacks=[es],
                )
            train_crps = model.evaluate([tfc_mean_train, tfc_std_train], 
                                        tobs_train, verbose=0)[0]
            valid_crps = model.evaluate([tfc_mean_test, tfc_std_test], 
                                        tobs_test, verbose=0)[0]
            if verbose == 1:
                print(train_crps, valid_crps)
        else:
            # For the more general network combine the inputs
            in_train = np.column_stack([tfc_mean_train, tfc_std_train])
            in_test = np.column_stack([tfc_mean_test, tfc_std_test])
            es = EarlyStopping(monitor='loss', min_delta=early_stopping_delta, 
                               patience=2)
            batch_size=in_train.shape[0]
            model.fit(
                in_train, tobs_train,
                epochs=epochs_max,
                batch_size=batch_size,
                verbose=verbose,
                callbacks=[es],
                )
            train_crps = model.evaluate(in_train, tobs_train, verbose=0)[0]
            valid_crps = model.evaluate(in_test, tobs_test, verbose=0)[0]

        # Write output
        train_crps_list.append(train_crps)
        valid_crps_list.append(valid_crps)

    # Stop timer 
    time_stop = timeit.default_timer()
    print('Time: %.2f s' % (time_stop - time_start))

    return train_crps_list, valid_crps_list


def load_data(data_dir, aux_dict=None):
    """Load data from NetCDF files.
    
    Ensemble mean and standard deviations are appended to feature list.
    Except for geo variables which do not have an ensemble dimension.

    Params:
        data_dir: base directory where NetCDF files are stored
        aux_dict: Dictionary with name of auxiliary file and 
                  list of variables. If None, only temperature.
    Returns:
        target: t2m_obs [date, station]
        feature_list: [feature, date, station]
        dates: [dates] List of datetime objects
    """
    aux_dir = data_dir + 'auxiliary/interpolated_to_stations/'
    
    fl = []   # Here we will store all the features
    # Load Temperature data
    rg = Dataset(data_dir + 'data_interpolated_00UTC.nc')
    target = rg.variables['t2m_obs'][:]
    ntime = target.shape[0]
    dates = num2date(rg.variables['time'][:],
                     units=rg.variables['time'].units)
    station_id = rg.variables['station_id'][:]
    fl.append(np.mean(rg.variables['t2m_fc'][:], axis=1))
    fl.append(np.std(rg.variables['t2m_fc'][:], axis=1, ddof=1))
    rg.close()
    
    if aux_dict is not None:
        for aux_fn, var_list in aux_dict.items():
            rg = Dataset(aux_dir + aux_fn)
            for var in var_list:
                data = rg.variables[var][:]
                if 'geo' in aux_fn:   
                    # Should probably look at dimensions
                    fl.append(np.array([data] * ntime))
                else:
                    fl.append(np.mean(data, axis=1))
                    fl.append(np.std(data, axis=1, ddof=1))
            rg.close()
    
    return target, np.array(fl), dates, station_id


def split_and_scale(target, features, dates, station_id, train_date_idx_start, 
                    train_date_idx_stop, test_date_idx_start, 
                    test_date_idx_stop): 
    """Splits the dataset into train and test set. Then the features
    are scaled by dividing by the training set max.

    Params:
        target: [date, station]
        features: [feature, date, station]
        train_date_idx_start: date id where training set starts (inc)
        train_date_idx_stop: date id where training set stops (excl)
        test_date_idx_stop: date id where test set starts (inc)
        test_date_idx_stop: date id where test set stops (excl)

    Returns:
        features_train: [feature, instance]
        target_train: [instance]
        features_test: [feature, instance]
        target_test: [instance]
        id_array_train: [instance]  Containing continuous IDs for embedding
        id_array_test: [instance]
    
    """

    # Split data set
    print('Train set contains %i days' % 
          (train_date_idx_stop - train_date_idx_start))
    print('Test set contains %i days' % 
          (test_date_idx_stop - test_date_idx_start))
    features_train = features[:, train_date_idx_start:train_date_idx_stop]
    target_train = target[train_date_idx_start:train_date_idx_stop]
    features_test = features[:, test_date_idx_stop:test_date_idx_stop]
    target_test = target[test_date_idx_stop:test_date_idx_stop]

    # Ravel arrays
    features_train = np.reshape(features_train, (features_train.shape[0], -1))
    target_train = np.reshape(target_train, (-1))
    features_test = np.reshape(features_test, (features_train.shape[0], -1))
    target_test = np.reshape(target_test, (-1))

    # Remove nans
    train_mask = np.isfinite(target_train.data)
    features_train = features_train[:, train_mask]
    target_train = target_train.data[train_mask]
    test_mask = np.isfinite(target_test.data)
    features_test = features_test[:, test_mask]
    target_test = target_test.data[test_mask]

    # Swap axes
    features_train = np.rollaxis(features_train, 1, 0)
    features_test = np.rollaxis(features_test, 1, 0)

    # Scale features
    features_max = np.max(features_train, axis=0)
    target_max = np.max(target_train)
    features_train /= features_max
    features_test /= features_max
    # target_train /= target_max   # No scaling of the outputs!
    # target_test /= target_max
    
    # Create continuous id array
    s = features.shape
    id_array = np.array([np.arange(s[-1])] * s[1])
    id_array_train = id_array[train_date_idx_start:train_date_idx_stop]
    id_array_test = id_array[test_date_idx_stop:test_date_idx_stop]
    id_array_train = np.reshape(id_array_train, (-1))
    id_array_test = np.reshape(id_array_test, (-1))
    id_array_train = id_array_train[train_mask]
    id_array_test = id_array_test[test_mask]

    # Create actual station id array
    station_array = np.array([list(station_id)] * s[1])
    station_array_train = station_array[train_date_idx_start:train_date_idx_stop]
    station_array_test = station_array[test_date_idx_stop:test_date_idx_stop]
    station_array_train = np.reshape(station_array_train, (-1))
    station_array_test = np.reshape(station_array_test, (-1))
    station_array_train = station_array_train[train_mask]
    station_array_test = station_array_test[test_mask]

    # Creat date array
    date_str = [datetime.strftime(dt, '%Y-%m-%d') for dt in list(dates)]
    date_array = np.array([date_str] * s[2])
    date_array = np.rollaxis(date_array, 1, 0)
    date_array_train = date_array[train_date_idx_start:train_date_idx_stop]
    date_array_test = date_array[test_date_idx_stop:test_date_idx_stop]
    date_array_train = np.reshape(date_array_train, (-1))
    date_array_test = np.reshape(date_array_test, (-1))
    date_array_train = date_array_train[train_mask]
    date_array_test = date_array_test[test_mask]
    
    train_set = Data(target_train, features_train, id_array_train,
                        station_array_train, date_array_train)
    test_set = Data(target_test, features_test, id_array_test,
                        station_array_test, date_array_test)
    return train_set, test_set


class Data(object):
    """Class for storing data
    """
    def __init__(self, target, features, cont_id, station_id, date_strs):
        self.target = target
        self.features = features
        self.cont_id = cont_id
        self.station_id = station_id
        self.date_strs = date_strs


