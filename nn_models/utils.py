"""
Utility functions for PPNN

Author: Stephan Rasp
"""

from scipy.stats import norm
import numpy as np
from netCDF4 import num2date, Dataset


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






