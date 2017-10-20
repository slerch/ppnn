"""
Utility functions for PPNN.
The functions are still a little confusingly named and 
structured at the moment.

Author: Stephan Rasp
"""
import os
print('Anaconda environment:', os.environ['CONDA_DEFAULT_ENV'])

from scipy.stats import norm
import numpy as np
from netCDF4 import num2date, Dataset
from emos_network_theano import EMOS_Network
import timeit
from keras.callbacks import EarlyStopping
from datetime import datetime
from tqdm import tqdm
import pandas as pd
import pdb

# Basic setup
date_format = '%Y-%m-%d'

# Data loading functions
def get_train_test_sets(data_dir=None, train_dates=None, test_dates=None,
                        predict_date=None, fclt=None, window_size=None,
                        preloaded_data=None, aux_dict=None,
                        verbose=1, seq_len=None, fill_value=None):
    """Load data and return train and test set objects.
    
    Parameters:
        data_dir: base directory where NetCDF files are stored
        train_dates_idxs: list with start (inc) and stop (excl) date str
                          yyyy-mm-dd
        test_dates_idxs: list with start (inc) and stop (excl) date str
        preloaded_data: list with [target, feature_list, dates, station_id, 
                        feature_names]
        aux_dict: Dictionary with name of auxiliary file and 
                  list of variables. If None, only temperature.
        seq_len: If given, return sequence data for RNNs
        fill_value: If given, convert missing data to fill_value
    """

    # Load raw data from netcdf files
    if preloaded_data is None:
        raw_data = load_raw_data(data_dir, aux_dict)
    else:
        raw_data = preloaded_data

    # Get date indices
    if train_dates is not None:
        train_dates_idxs = [return_date_idx(raw_data[2], train_dates[0]),
                            return_date_idx(raw_data[2], train_dates[1])]
        test_dates_idxs = [return_date_idx(raw_data[2], test_dates[0]),
                           return_date_idx(raw_data[2], test_dates[1])]
    else:
        date_idx = return_date_idx(raw_data[2], predict_date)
        fclt_didx = int(fclt / 24)
        train_dates_idxs = [date_idx - fclt_didx - window_size, 
                            date_idx - fclt_didx]
        test_dates_idxs = [date_idx, date_idx + 1]

    # Split into test and train set and scale features
    train_set, test_set = split_and_scale(raw_data, train_dates_idxs, 
                                          test_dates_idxs, verbose,
                                          seq_len, fill_value)

    return train_set, test_set


def load_raw_data(data_dir, aux_dict=None):
    """Load raw data from NetCDF files.
    
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
    feature_names = ['t2m_fc_mean', 't2m_fc_std']
    
    if aux_dict is not None:
        for aux_fn, var_list in aux_dict.items():
            rg = Dataset(aux_dir + aux_fn)
            for var in var_list:
                data = rg.variables[var][:]
                if 'geo' in aux_fn:   
                    # Should probably look at dimensions
                    fl.append(np.array([data] * ntime))
                    feature_names.extend(var)
                else:
                    fl.append(np.mean(data, axis=1))
                    fl.append(np.std(data, axis=1, ddof=1))
                    feature_names.extend([var + '_mean', var + '_std'])
            rg.close()
    
    return (target.data, np.array(fl, dtype='float32'), dates, station_id, 
            feature_names)


class DataContainer(object):
    """Class for storing data
    """
    def __init__(self, targets, features, cont_ids, station_ids, date_strs,
                 feature_names, sample_weights=None):
        self.targets = targets
        self.features = features
        self.cont_ids = cont_ids
        self.station_ids = station_ids
        self.date_strs = date_strs
        self.feature_names = feature_names
        self.sample_weights = sample_weights


def split_and_scale(raw_data, train_dates_idxs, test_dates_idxs, verbose=1,
                    seq_len=None, fill_value=None):
    """
    """

    # Unpack raw_data
    targets, features, dates, station_id, feature_names = raw_data

    data_sets = []
    for set_name, dates_idxs in zip(['train', 'test'], 
                                    [train_dates_idxs, test_dates_idxs]):

        # Split data set: 
        if verbose == 1:
            print('%s set contains %i days' % 
                  (set_name, dates_idxs[1] - dates_idxs[0]))

        if seq_len is None:
            t = targets[dates_idxs[0]:dates_idxs[1]] # [date, station]
            f = features[:, dates_idxs[0]:dates_idxs[1]] # [feature, date, station]

            # Ravel arrays, combine dates and stations --> instances
            t = np.reshape(t, (-1)) # [instances]
            f = np.reshape(f, (f.shape[0], -1)) # [features, instances]

            # Swap feature axes
            f = np.rollaxis(f, 1, 0) # [instances, features]

            # Get nan mast from target
            nan_mask = np.isfinite(t)
        else:
            t = targets[dates_idxs[0]-seq_len+1:dates_idxs[1]] # [date, station]
            f = features[:, dates_idxs[0]-seq_len+1:dates_idxs[1]] # [feature, date, station]

            # Stack time steps for sequences
            # [time_step, feature, day, station]
            t = np.stack([t[i:-(seq_len-i-1) or None] for i in range(seq_len)])
            # [time_step, day, station]
            f = np.stack([f[:, i:-(seq_len-i-1) or None] for i in range(seq_len)])

            # Ravel arrays [seq, feature, instance]
            t = np.reshape(t, (seq_len, -1))
            f = np.reshape(f, (seq_len, f.shape[1], -1))

            # Roll arrays[sample, time step, feature]
            t = np.rollaxis(t, 1, 0)
            f = np.rollaxis(f, 2, 0)
            t = np.atleast_3d(t)

            # Get nan mask from last entry of target
            nan_mask = np.isfinite(t[:, -1, 0])
        
        # Apply NaN mask
        f = f[nan_mask]
        t = t[nan_mask]

        # Scale features
        if set_name == 'train': # Get maximas
            if seq_len is None:
                features_max = np.max(f, axis=0)
            else:
                features_max = np.max(f, axis=(0, 1))
        f /= features_max

        # Replace NaNs with fill value is requested
        if fill_value is not None:
            assert seq_len is not None, 'fill value only implemented for sequences.'
            weights = np.array(np.isfinite(t[:, :, 0]), dtype=np.float32)
            t[np.isnan(t)] = fill_value
        else:
            weights = None
        
        # Get additional data
        cont_ids = get_cont_ids(features, nan_mask, dates_idxs)
        station_ids = get_station_ids(features, station_id, nan_mask, dates_idxs)
        date_strs = get_date_strs(features, dates, nan_mask, dates_idxs)

        # Put in data container
        data_sets.append(DataContainer(t, f, cont_ids, station_ids, 
                                       date_strs, feature_names, weights))

    return data_sets


def get_cont_ids(features, nan_mask, dates_idxs):
    s = features.shape
    cont_ids = np.array([np.arange(s[-1])] * s[1])
    cont_ids = cont_ids[dates_idxs[0]:dates_idxs[1]]
    cont_ids = np.reshape(cont_ids, (-1))
    return cont_ids[nan_mask]


def get_station_ids(features, station_id, nan_mask, dates_idxs):
    s = features.shape
    station_ids = np.array([list(station_id)] * s[1])
    station_ids = station_ids[dates_idxs[0]:dates_idxs[1]]
    station_ids = np.reshape(station_ids, (-1))
    return station_ids[nan_mask]


def get_date_strs(features, dates, nan_mask, dates_idxs):
    s = features.shape
    date_strs = [datetime.strftime(dt, date_format) for dt in list(dates)]
    date_strs = np.array([list(date_strs)] * s[2])
    date_strs = np.rollaxis(date_strs, 1, 0)
    date_strs = date_strs[dates_idxs[0]:dates_idxs[1]]
    date_strs = np.reshape(date_strs, (-1))
    return date_strs[nan_mask]


# Helper functions
def return_date_idx(dates, date_str=None, y=None, m=None, d=None):
    if date_str is not None:
        dt = datetime.strptime(date_str, date_format)
    else:
        dt = datetime(y, m, d, 0, 0)
    return np.where(dates == dt)[0][0]


def crps_normal(mu, sigma, y):
    """
    Compute CRPS for a Gaussian distribution. 
    """
    loc = (y - mu) / sigma
    crps = sigma * (loc * (2 * norm.cdf(loc) - 1) + 
                    2 * norm.pdf(loc) - 1. / np.sqrt(np.pi))
    return crps


def loop_over_days(data_dir, model, date_str_start, date_str_stop, 
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
    date_str_start/stop : int
        Start and stop string for loop
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

    # Load data initially
    raw_data = load_raw_data(data_dir)

    # Get date indices
    dates = raw_data[2]
    date_idx_start = return_date_idx(dates, date_str_start)
    date_idx_stop = return_date_idx(dates, date_str_stop)
    date_str_list = [datetime.strftime(dt, date_format) for dt in 
                     list(dates[date_idx_start:date_idx_stop])]


    # Initialize lists to store dates, station_ids, means and stds
    date_list = []
    station_id_list =[]
    mean_list = []
    std_list = [] 

    # Start loop 
    for i, date_str in enumerate(tqdm(date_str_list)):

        # Get data slice
        train_set, test_set = get_train_test_sets(preloaded_data=raw_data, 
                                                  predict_date=date_str,
                                                  fclt=fclt, 
                                                  window_size=window_size,
                                                  verbose=verbose)

        # Write dates and station_ids in list
        date_list.extend(list(test_set.date_strs))
        station_id_list.extend(list(test_set.station_ids))

        # Reinitialize model if requested
        # only works for theano model
        if reinit_model:
            assert model_type == 'EMOS_Network_theano', \
                'Reinitialize does not work with keras'
            model = EMOS_Network()

        # Fit model
        if model_type == 'EMOS_Network_theano':
            # Split mean and std
            train_mean = train_set.features[:, 0]
            train_std = train_set.features[:, 1]
            test_mean = test_set.features[:, 0]
            test_std = test_set.features[:, 1]
            # Fit model
            train_crps, valid_crps = model.fit(
                train_mean, train_std, train_set.targets, 
                epochs_max, 
                (test_mean, test_std, test_set.targets), 
                lr=lr, 
                early_stopping_delta=early_stopping_delta,
                verbose=verbose,
                )
            # Get prediction
            p = model.predict(test_mean, test_std, test_set.targets)
            p_means, p_stds = p[0], p[1]

        elif model_type == 'EMOS_Network_keras':
            # Split mean and std
            train_mean = train_set.features[:, 0]
            train_std = train_set.features[:, 1]
            test_mean = test_set.features[:, 0]
            test_std = test_set.features[:, 1]
            # Setup
            es = EarlyStopping(monitor='loss', min_delta=early_stopping_delta, 
                               patience=2)
            batch_size=train_mean.shape[0]
            # Fit model
            model.fit(
                [train_mean, train_std], train_set.targets,
                epochs=epochs_max,
                batch_size=batch_size,
                verbose=verbose,
                callbacks=[es],
                )
            train_crps = model.evaluate([train_mean, train_std], 
                                        train_set.targets, verbose=0)[0]
            valid_crps = model.evaluate([test_mean, test_std], 
                                        test_set.targets, verbose=0)[0]
            if verbose == 1:
                print(train_crps, valid_crps)

            # Get predictions
            p = model.predict([test_mean, test_std])
        else:
            # For the more general network combine the inputs
            es = EarlyStopping(monitor='loss', min_delta=early_stopping_delta, 
                               patience=2)
            batch_size=train_set.features.shape[0]
            model.fit(
                train_set.features, train_set.targets,
                epochs=epochs_max,
                batch_size=batch_size,
                verbose=verbose,
                callbacks=[es],
                )
            train_crps = model.evaluate(train_set.features, train_set.targets, verbose=0)[0]
            valid_crps = model.evaluate(test_set.features, test_set.targets, verbose=0)[0]

            # Get predictions
            p = model.predict(test_set.features)

        # Write output
        train_crps_list.append(train_crps)
        valid_crps_list.append(valid_crps)

        # Store predictions
        p_means, p_stds = p[:, 0], p[:, 1]
        mean_list.extend(list(p_means))
        std_list.extend(list(p_stds))

    # Create pandas dataframe
    results_df = create_results_df(date_list, station_id_list, mean_list, 
                                   std_list)

    return train_crps_list, valid_crps_list, results_df


def create_results_df(dates, station_ids, means, stds):
    """
    """
    df = pd.DataFrame({
        'date': dates,
        'station_id': station_ids,
        'mean': means,
        'std': stds,
        })
    return df


# def split_and_scale(target, features, dates, station_id, train_date_idx_start, 
#                     train_date_idx_stop, test_date_idx_start, 
#                     test_date_idx_stop): 
#     """Splits the dataset into train and test set. Then the features
#     are scaled by dividing by the training set max.

#     Params:
#         target: [date, station]
#         features: [feature, date, station]
#         train_date_idx_start: date id where training set starts (inc)
#         train_date_idx_stop: date id where training set stops (excl)
#         test_date_idx_stop: date id where test set starts (inc)
#         test_date_idx_stop: date id where test set stops (excl)

#     Returns:
#         features_train: [feature, instance]
#         target_train: [instance]
#         features_test: [feature, instance]
#         target_test: [instance]
#         id_array_train: [instance]  Containing continuous IDs for embedding
#         id_array_test: [instance]
    
#     """

#     # Split data set
#     print('Train set contains %i days' % 
#           (train_date_idx_stop - train_date_idx_start))
#     print('Test set contains %i days' % 
#           (test_date_idx_stop - test_date_idx_start))
#     features_train = features[:, train_date_idx_start:train_date_idx_stop]
#     target_train = target[train_date_idx_start:train_date_idx_stop]
#     features_test = features[:, test_date_idx_stop:test_date_idx_stop]
#     target_test = target[test_date_idx_stop:test_date_idx_stop]

#     # Ravel arrays
#     features_train = np.reshape(features_train, (features_train.shape[0], -1))
#     target_train = np.reshape(target_train, (-1))
#     features_test = np.reshape(features_test, (features_train.shape[0], -1))
#     target_test = np.reshape(target_test, (-1))

#     # Remove nans
#     train_mask = np.isfinite(target_train.data)
#     features_train = features_train[:, train_mask]
#     target_train = target_train.data[train_mask]
#     test_mask = np.isfinite(target_test.data)
#     features_test = features_test[:, test_mask]
#     target_test = target_test.data[test_mask]

#     # Swap axes
#     features_train = np.rollaxis(features_train, 1, 0)
#     features_test = np.rollaxis(features_test, 1, 0)

#     # Scale features
#     features_max = np.max(features_train, axis=0)
#     target_max = np.max(target_train)
#     features_train /= features_max
#     features_test /= features_max
#     # target_train /= target_max   # No scaling of the outputs!
#     # target_test /= target_max
    
#     # Create continuous id array
#     s = features.shape
#     id_array = np.array([np.arange(s[-1])] * s[1])
#     id_array_train = id_array[train_date_idx_start:train_date_idx_stop]
#     id_array_test = id_array[test_date_idx_stop:test_date_idx_stop]
#     id_array_train = np.reshape(id_array_train, (-1))
#     id_array_test = np.reshape(id_array_test, (-1))
#     id_array_train = id_array_train[train_mask]
#     id_array_test = id_array_test[test_mask]

#     # Create actual station id array
#     station_array = np.array([list(station_id)] * s[1])
#     station_array_train = station_array[train_date_idx_start:train_date_idx_stop]
#     station_array_test = station_array[test_date_idx_stop:test_date_idx_stop]
#     station_array_train = np.reshape(station_array_train, (-1))
#     station_array_test = np.reshape(station_array_test, (-1))
#     station_array_train = station_array_train[train_mask]
#     station_array_test = station_array_test[test_mask]

#     # Creat date array
#     date_str = [datetime.strftime(dt, '%Y-%m-%d') for dt in list(dates)]
#     date_array = np.array([date_str] * s[2])
#     date_array = np.rollaxis(date_array, 1, 0)
#     date_array_train = date_array[train_date_idx_start:train_date_idx_stop]
#     date_array_test = date_array[test_date_idx_stop:test_date_idx_stop]
#     date_array_train = np.reshape(date_array_train, (-1))
#     date_array_test = np.reshape(date_array_test, (-1))
#     date_array_train = date_array_train[train_mask]
#     date_array_test = date_array_test[test_mask]
    
#     train_set = Data(target_train, features_train, id_array_train,
#                         station_array_train, date_array_train)
#     test_set = Data(target_test, features_test, id_array_test,
#                         station_array_test, date_array_test)
#     return train_set, test_set


# class Data(object):
#     """Class for storing data
#     """
#     def __init__(self, target, features, cont_id, station_id, date_strs):
#         self.target = target
#         self.features = features
#         self.cont_id = cont_id
#         self.station_id = station_id
#         self.date_strs = date_strs


