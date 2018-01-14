"""
Utility functions for PPNN.
The functions are still a little confusingly named and 
structured at the moment.

Author: Stephan Rasp
"""
import os
import platform
import copy
from scipy.stats import norm
import numpy as np
np.random.seed(42)
from netCDF4 import num2date, Dataset
from emos_network_theano import EMOS_Network
import timeit
from keras.callbacks import EarlyStopping
from datetime import datetime
from tqdm import tqdm
import pandas as pd
from collections import OrderedDict
import pdb
import pickle
import matplotlib.pyplot as plt

# Basic setup
print('Anaconda environment:', os.environ['CONDA_DEFAULT_ENV'])
print(platform.system(), platform.release())
date_format = '%Y-%m-%d'


# Define aux variable dictionary
aux_dict = OrderedDict()
aux_dict['data_aux_geo_interpolated.nc'] = ['orog', 
                                            'station_alt', 
                                            'station_lat', 
                                            'station_lon']
aux_dict['data_aux_pl500_interpolated_00UTC.nc'] = ['u_pl500_fc',
                                                    'v_pl500_fc',
                                                    'gh_pl500_fc']
aux_dict['data_aux_pl850_interpolated_00UTC.nc'] = ['u_pl850_fc',
                                                    'v_pl850_fc',
                                                    'q_pl850_fc']
aux_dict['data_aux_surface_interpolated_00UTC.nc'] = ['cape_fc',
                                                      'sp_fc',
                                                      'tcc_fc']
aux_dict['data_aux_surface_more_interpolated_part1_00UTC.nc']  = [
    'sshf_fc', 'slhf_fc', 'u10_fc','v10_fc'
]
aux_dict['data_aux_surface_more_interpolated_part2_00UTC.nc']  = [
    'ssr_fc', 'str_fc', 'd2m_fc','sm_fc'
]


# Data loading functions
def get_train_test_sets(data_dir=None, train_dates=None, test_dates=None,
                        predict_date=None, fclt=48, window_size=None,
                        preloaded_data=None, aux_dict=None,
                        verbose=1, seq_len=None, fill_value=None,
                        valid_size=None, full_ensemble_t=False,
                        add_current_error=False, current_error_len=1):
    """Load data and return train and test set objects.
    
    Args:
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
        valid_size: If given, returns a third set containing a given fraction of
                    the train set.
        full_ensemble_t: Return all 50 ensemble members for temperature
        add_current_error: If True, add current mean forecast, obs and error
        fclt: Forecast lead time in hours. Used if add_current_error is True.
              Default = 48h
    Returns:
        train_set, test_set: Objects containing train and test data
    """

    # Load raw data from netcdf files
    if preloaded_data is None:
        raw_data = load_raw_data(data_dir, aux_dict, full_ensemble_t)
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
                                          seq_len, fill_value, 
                                          full_ensemble_t, add_current_error,
                                          fclt, current_error_len)
    
    # Split train set if requested
    if valid_size is not None:
        arrays = [train_set.targets, train_set.features, train_set.cont_ids, 
                  train_set.station_ids, train_set.date_strs, 
                  train_set.sample_weights]
        train_arrays, valid_arrays = ([], [])

        random_idxs = np.arange(arrays[0].shape[0])
        np.random.shuffle(random_idxs)
        split_idx = int(random_idxs.shape[0] * 0.2)
        train_idxs = random_idxs[split_idx:]
        valid_idxs = random_idxs[:split_idx]

        valid_set = copy.deepcopy(train_set)
        for s, idxs in zip([train_set, valid_set], [train_idxs, valid_idxs]):
            s.targets = arrays[0][idxs]
            s.features = arrays[1][idxs]
            s.cont_ids = arrays[2][idxs]
            s.station_ids = arrays[3][idxs]
            s.date_strs = arrays[4][idxs]
            s.sample_weights = arrays[5][idxs]
        return train_set, test_set, valid_set
    else:
        return train_set, test_set


def load_raw_data(data_dir, aux_dict=None, full_ensemble_t=False):
    """Load raw data from NetCDF files.
    
    Ensemble mean and standard deviations are appended to feature list.
    Except for geo variables which do not have an ensemble dimension.

    Params:
        data_dir: base directory where NetCDF files are stored
        aux_dict: Dictionary with name of auxiliary file and 
                  list of variables. If None, only temperature.
        full_ensemble_t: Return all 50 ensemble members for temperature
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
    if full_ensemble_t:
        feature_names = []
        for i in range(rg.variables['t2m_fc'].shape[1]):
            fl.append(rg.variables['t2m_fc'][:, i, :])
            feature_names.append('t2m_fc_ens%i' % (i+1))
    else:
        fl.append(np.mean(rg.variables['t2m_fc'][:], axis=1))
        fl.append(np.std(rg.variables['t2m_fc'][:], axis=1, ddof=1))
        feature_names = ['t2m_fc_mean', 't2m_fc_std']
    rg.close()
    
    if aux_dict is not None:
        for aux_fn, var_list in aux_dict.items():
            rg = Dataset(aux_dir + aux_fn)
            for var in var_list:
                data = rg.variables[var][:]
                if 'geo' in aux_fn:   
                    # Should probably look at dimensions
                    fl.append(np.array([data] * ntime))
                    feature_names.extend([var])
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
                 feature_names, sample_weights=None, scale_factors=None):
        self.targets = targets
        self.features = features
        self.cont_ids = cont_ids
        self.station_ids = station_ids
        self.date_strs = date_strs
        self.feature_names = feature_names
        self.sample_weights = sample_weights
        self.scale_factors = scale_factors


def split_and_scale(raw_data, train_dates_idxs, test_dates_idxs, verbose=1,
                    seq_len=None, fill_value=None, full_ensemble_t=False,
                    add_current_error=False, fclt=48, current_error_len=1):
    """
    """

    # Unpack raw_data
    targets, features, dates, station_id, feature_names = raw_data

    if add_current_error:
        feature_names.extend(['curr_t2m_fc_obs',
                              'curr_err'])
        if current_error_len > 1:
            for i in range(1, current_error_len, 1):
                feature_names.extend([
                    'curr_t2m_fc_obs_m%i' % i,
                    'curr_err_m%i' % i
                ])
        assert full_ensemble_t is False, 'Current error not compatible with full ensemble.'

    data_sets = []
    for set_name, dates_idxs in zip(['train', 'test'], 
                                    [train_dates_idxs, test_dates_idxs]):

        # Split data set: 
        if verbose == 1:
            print('%s set contains %i days' % 
                  (set_name, dates_idxs[1] - dates_idxs[0]))

        if seq_len is None:
            #pdb.set_trace()
            t = targets[dates_idxs[0]:dates_idxs[1]] # [date, station]
            f = features[:, dates_idxs[0]:dates_idxs[1]] # [feature, date, station]

            if add_current_error:
                didx = int(fclt / 24)
                new_f_list = []
                for i in range(current_error_len):
                    d = didx + i
                    curr_obs = targets[dates_idxs[0]-d:dates_idxs[1]-d].copy()
                    curr_fc = features[0, dates_idxs[0] - d:dates_idxs[1] - d]
                    # Replace missing observations with forecast values
                    # [date_shifted, station]
                    curr_obs[np.isnan(curr_obs)] = curr_fc[np.isnan(curr_obs)]
                    curr_err = curr_obs - curr_fc
                    new_f_list.extend([curr_obs, curr_err])
    
                new_f = np.stack(new_f_list, axis=0)
                # [new features, date_shifted, station]
                f = np.concatenate((f, new_f), axis=0)

            # Ravel arrays, combine dates and stations --> instances
            t = np.reshape(t, (-1)) # [instances]
            f = np.reshape(f, (f.shape[0], -1)) # [features, instances]

            # Swap feature axes
            f = np.rollaxis(f, 1, 0) # [instances, features]

            # Get nan mask from target
            nan_mask = np.isfinite(t)
        else:
            assert add_current_error is False, 'Current error not compatible with sequence'
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
                features_max = np.nanmax(f, axis=0)
            else:
                features_max = np.nanmax(f, axis=(0, 1))
        if full_ensemble_t:
            # Scale all temeperature members with same max
            n_ens = 50  # ATTENTION: hard-coded
            features_max[:n_ens] = np.max(features_max[:n_ens])
        f /= features_max

        # Replace NaNs with fill value is requested
        if fill_value is not None:
            assert seq_len is not None, 'fill value only implemented for sequences.'
            weights = np.array(np.isfinite(t[:, :, 0]), dtype=np.float32)
            t[np.isnan(t)] = fill_value
        else:
            weights = None
        
        # Get additional data
        cont_ids = get_cont_ids(features, nan_mask, dates_idxs, seq_len)
        station_ids = get_station_ids(features, station_id, nan_mask, dates_idxs)
        date_strs = get_date_strs(features, dates, nan_mask, dates_idxs)

        # Put in data container
        data_sets.append(DataContainer(t, f, cont_ids, station_ids, 
                                       date_strs, feature_names, weights,
                                       features_max))

    return data_sets


def get_cont_ids(features, nan_mask, dates_idxs, seq_len=None):
    s = features.shape
    cont_ids = np.array([np.arange(s[-1])] * s[1])
    cont_ids = cont_ids[dates_idxs[0]:dates_idxs[1]]
    cont_ids = np.reshape(cont_ids, (-1))
    if seq_len is not None:
        cont_ids = np.rollaxis(np.array([cont_ids] * seq_len), 1, 0)
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


def maybe_correct_cat_crps(preds, targets, bin_edges):
    """CRPS for categorical predictions. I think this is correct now.

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
    adj_cum_probs = np.concatenate((np.zeros((cum_probs.shape[0], 1)),
                                    cum_probs), axis=1)
    # Compute squared area for each bin
    sq_list = []
    for i in range(cum_bin_obs.shape[1]):
        x_l = np.abs(cum_bin_obs[:, i] - adj_cum_probs[:, i])
        x_r = np.abs(cum_bin_obs[:, i] - adj_cum_probs[:, i + 1])
        sq = 1./3. * (x_l ** 2 + x_l * x_r + x_r ** 2)
        sq_list.append(sq)

    # Compute CRPS
    crps = np.sum(np.array(sq_list).T * np.diff(ins_bin_edges, axis=1), axis=1)
    return np.mean(crps)


# Experiment running functions
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
            # Convert to keras format
            p = np.stack([p[0], p[1]], axis=1)

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
                                        train_set.targets, verbose=0)
            valid_crps = model.evaluate([test_mean, test_std], 
                                        test_set.targets, verbose=0)
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
            train_crps = model.evaluate(train_set.features, train_set.targets,
                                        verbose=0)
            valid_crps = model.evaluate(test_set.features, test_set.targets,
                                        verbose=0)

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


def create_results_df(dates, station_ids, means, stds, train_time=None,
                      params=None):
    """
    """
    df = pd.DataFrame({
        'date': dates,
        'station_id': station_ids,
        'mean': means,
        'std': stds,
        'train_time': train_time,
        'params': params
        })
    return df


def save_pickle(data_dir, fn, train_dates=['2015-01-01', '2016-01-01'],
                add_current_error=False, current_error_len=1,
                test_dates=['2016-01-01', '2017-01-01']):
    """Load and pickle dataset"""
    sets = get_train_test_sets(
        data_dir, train_dates, test_dates, aux_dict=aux_dict,
        add_current_error=add_current_error, current_error_len=current_error_len
    )
    with open(data_dir + fn, 'wb') as f:
        pickle.dump(sets, f)


def plot_fc(data_set, idx, distr='pdf', preds=None):
    fc = data_set.features[idx, :2] * data_set.scale_factors
    obs = data_set.targets[idx]

    x = np.linspace(fc[0] - 5 * fc[1], fc[0] + 5 * fc[1], 100)
    if distr == 'pdf':
        y = norm.pdf(x, fc[0], fc[1])
    elif distr == 'cdf':
        y = norm.cdf(x, fc[0], fc[1])
    else:
        raise Exception
    plt.plot(x, y, label='raw ensemble')
    plt.axvline(obs, color='red', label='obs')
    if preds is not None:
        p = preds[idx]
        if distr == 'pdf':
            y = norm.pdf(x, p[0], p[1])
        elif distr == 'cdf':
            y = norm.cdf(x, p[0], p[1])
        plt.plot(x, y, label='prediction')
    plt.xlabel('Temperature [C]')
    plt.legend()
    plt.show()