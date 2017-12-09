"""Runs network experiment

"""

from configargparse import ArgParser
import numpy as np
from utils import get_train_test_sets, create_results_df
from keras_models import build_fc_model, build_hidden_model, build_emb_model
import keras
from keras.callbacks import EarlyStopping
from losses import crps_cost_function
from aux_dict import aux_dict
import pdb
import pickle


def main(inargs):
    """Main program to run network experiment

    Args:
        inargs:

    Returns:

    """
    # Load data
    if verbose > 0: print('Loading train and test data.')
    if inargs.pickled_sets is None:
        var_dict = aux_dict if inargs.use_aux else None
        train_set, test_set = get_train_test_sets(
            inargs.data_dir,
            inargs.train_dates,
            inargs.test_dates,
            aux_dict=var_dict,
            add_current_error=inargs.add_current_error,
            current_error_len=inargs.current_error_len,
        )
    else:
        with open(inargs.pickled_sets, 'rb') as f:
            train_set, test_set = pickle.load(f)

    # Build keras model
    n_features = train_set.features.shape[1]
    n_outputs = 2
    if inargs.model == 'fc':
        model = build_fc_model(
            n_features,
            n_outputs,
            compile=True,
            lr=inargs.lr,
        )
    elif inargs.model == 'hidden':
        model = build_hidden_model(
            n_features,
            n_outputs,
            hidden_nodes=inargs.hidden_nodes,
            compile=True,
            lr=inargs.lr,
        )
    elif inargs.model == 'emb':
        max_id = int(np.max([train_set.cont_ids.max(),
                             test_set.cont_ids.max()]))
        model = build_emb_model(
            n_features,
            n_outputs,
            hidden_nodes=inargs.hidden_nodes,
            emb_size=inargs.emb_size,
            max_id=max_id,
            compile=True,
            lr=inargs.lr,
        )
    else:
        raise ValueError('Wrong model type.')
    if verbose > 0: print(model.summary())

    # Compile model
    if inargs.es_patience is not None:
        callbacks = [
            EarlyStopping(
                monitor='val_loss',
                patience=inargs.es_patience,
            )
        ]
    else: callbacks = []

    # Train model
    if inargs.model == 'emb':
        train_features = [train_set.features, train_set.cont_ids]
        test_features = [test_set.features, test_set.cont_ids]
    else:
        train_features = train_set.features
        test_features = test_set.features
    model.fit(
        train_features,
        train_set.targets,
        epochs=inargs.epochs,
        batch_size=inargs.batch_size,
        validation_split=inargs.validation_split,
        verbose=verbose,
        callbacks=callbacks,
    )

    # Test model
    print('Test score:', model.evaluate(
        test_features,
        test_set.targets,
        batch_size=4096,
        verbose=0,
    ))

    # Save predictions
    if inargs.save_preds:
        if verbose > 0: print('Save predictions')
        preds = model.predict(test_features, batch_size=4096)
        results_df = create_results_df(
            test_set.date_strs,
            test_set.station_ids,
            preds[:, 0],
            preds[:, 1]
        )
        results_df.to_csv(inargs.results_dir + inargs.exp_name + '.csv')


if __name__ == '__main__':

    p = ArgParser()

    # Config file
    p.add(
        '-c', '--config',
        is_config_file=True,
        help='Config file path.'
    )

    # Directories and experiment name
    p.add_argument(
        '--data_dir',
        type=str,
        default='/project/meteo/w2w/C7/ppnn_data/',
        help='Directory containing input data. '
             'Default: /project/meteo/w2w/C7/ppnn_data/',
    )
    p.add_argument(
        '--results_dir',
        type=str,
        default='../results/csv_files/',
        help='Directory to save results to. '
             'Default: ../results/',
    )
    p.add_argument(
        '--exp_name',
        type=str,
        default=None,
        help='Experiment name. Default: None',
    )

    # Train and test dates, data settings
    p.add_argument(
        '--train_dates',
        type=str,
        nargs='+',
        default=['2015-01-01', '2016-01-01'],
        help='Training start (inc) and stop (exc) dates in format '
             'yyyy-mm-dd. Default: 2015-01-01 2016-01-01',
    )
    p.add_argument(
        '--test_dates',
        type=str,
        nargs='+',
        default=['2016-01-01', '2017-01-01'],
        help='Test start (inc) and stop (exc) dates in format '
             'yyyy-mm-dd. Default: 2016-01-01 2017-01-01',
    )
    p.add_argument(
        '--use_aux',
        dest='use_aux',
        action='store_true',
        help='If given, use auxiliary variables.',
    )
    p.set_defaults(use_aux=False)
    p.add_argument(
        '--add_current_error',
        dest='add_current_error',
        action='store_true',
        help='If given, use current error.',
    )
    p.set_defaults(add_current_error=False)
    p.add_argument(
        '--current_error_len',
        type=int,
        default=1,
        help='Length of current error. Default: 1',
    )
    p.set_defaults(current_error=False)
    p.add_argument(
        '--pickled_sets',
        type=str,
        default=None,
        help='Load pickled datasets. Default: None',
    )

    # Network parameters
    p.add_argument(
        '--model',
        type=str,
        help='Type of the model: [fc, hidden, emb]',
    )
    p.add_argument(
        '--epochs',
        type=int,
        default=10,
        help='Number of epochs. Default: 10',
    )
    p.add_argument(
        '--batch_size',
        type=int,
        default=512,
        help='Training batch size. Default: 512',
    )
    p.add_argument(
        '--validation_split',
        type=float,
        default=0.,
        help='Fraction of training data for validation. Default: 0.',
    )
    p.add_argument(
        '--hidden_nodes',
        type=int,
        nargs='+',
        default=[],
        help='Int or list of hidden nodes in hidden and emb networks. '
             'Default: [] means no hidden layers',
    )
    p.add_argument(
        '--es_patience',
        type=int,
        default=None,
        help='Number of rounds to wait for validation improvement. '
             'Default: None means no early stopping.',
    )
    p.add_argument(
        '--emb_size',
        type=int,
        default=3,
        help='Number of latent features. Default: 3',
    )
    p.add_argument(
        '--lr',
        type=float,
        default=0.001,
        help='Learning rate. Default: 0.001',
    )

    # Other settings
    p.add_argument(
        '--verbose',
        type=int,
        default=0,
        help='Verbosity level: 0 or 1. Default: 0',
    )
    p.add_argument(
        '--save_preds',
        dest='save_preds',
        action='store_true',
        help='If given, save predictions in results_dir with exp_name',
    )
    p.set_defaults(save_preds=False)
    args = p.parse_args()

    # Set global verbosity level
    global verbose
    verbose = args.verbose

    main(args)