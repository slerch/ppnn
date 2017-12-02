"""Runs network experiment

"""

from argparse import ArgumentParser
from utils import get_train_test_sets, create_results_df
from keras_models import build_fc_model


def main(inargs):
    """Main program to run network experiment

    Args:
        inargs:

    Returns:

    """
    # Load data
    if verbose > 0: print('Loading train and test data.')
    train_set, test_set = get_train_test_sets(
        inargs.data_dir,
        inargs.train_dates,
        inargs.test_dates,
    )

    # Build keras model
    n_features = train_set.features.shape[1]
    n_outputs = 2
    model = build_fc_model(
        n_features,
        n_outputs,
        compile=True,
    )
    if verbose > 0: print(model.summary())

    # Train model
    model.fit(
        train_set.features,
        train_set.targets,
        epochs=inargs.epochs,
        batch_size=inargs.batch_size,
        validation_split=inargs.validation_split,
        verbose=verbose,
    )

    # Test model
    print('Test score:', model.evaluate(
        test_set.features,
        test_set.targets,
        batch_size=4096,
        verbose=0,
    ))

    # Save predictions
    if inargs.save_preds:
        if verbose > 0: print('Save predictions')
        preds = model.predict(test_set.features, batch_size=4096)
        results_df = create_results_df(
            test_set.date_strs,
            test_set.station_ids,
            preds[:, 0],
            preds[:, 1]
        )
        results_df.to_csv(inargs.results_dir + inargs.exp_name + '.csv')


if __name__ == '__main__':

    p = ArgumentParser()

    # Directories and experiment name
    p.add_argument(
        '--data_dir',
        type=str,
        default='/Volumes/STICK/data/ppnn_data/',
        help='Directory containing input data.'
             'Default: /Volumes/STICK/data/ppnn_data/',
    )
    p.add_argument(
        '--results_dir',
        type=str,
        default='../results/',
        help='Directory to save results to.'
             'Default: ../results/',
    )
    p.add_argument(
        '--exp_name',
        type=str,
        default=None,
        help='Experiment name. Default: None',
    )

    # Train and test dates
    p.add_argument(
        '--train_dates',
        type=str,
        nargs='+',
        default=['2015-01-01', '2016-01-01'],
        help='Training start (inc) and stop (exc) dates in format'
             'yyyy-mm-dd. Default: 2015-01-01 2016-01-01',
    )
    p.add_argument(
        '--test_dates',
        type=str,
        nargs='+',
        default=['2016-01-01', '2017-01-01'],
        help='Test start (inc) and stop (exc) dates in format'
             'yyyy-mm-dd. Default: 2016-01-01 2017-01-01',
    )

    # Network parameters
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
        help='Verbosity level: 0 or 1. Default: 0',
    )
    p.set_defaults(save_preds=False)
    args = p.parse_args()

    # Set global verbosity level
    global verbose
    verbose = args.verbose

    main(args)