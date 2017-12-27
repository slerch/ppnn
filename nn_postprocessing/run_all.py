"""
Script to run all experiments defined by config files.
"""
from glob import glob
import os

if __name__ == '__main__':
    cwd = os.getcwd()
    if cwd == '/Users/stephanrasp/repositories/ppnn/nn_postprocessing':
        data_dir = '/Volumes/STICK/data/ppnn_data/'
    elif cwd == '/home/s/S.Rasp/repositories/ppnn/nn_postprocessing':
        data_dir = '/project/meteo/w2w/C7/ppnn_data/'
    else:
        raise Exception('Working directory not recognized.')
    config_files = sorted(glob('./config/*.yml'))
    for c in config_files:
        print('Running experiment:', c)
        os.system('python ./run_experiment.py -c %s --data_dir %s --save_preds'
                  % (c, data_dir))
