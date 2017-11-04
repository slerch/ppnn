# Neural Net Post-processing

All the neural net experiments are run in Jupyter notebooks which use functions defined in other files in this directory. 

| File | Description |
| ---- | ----------- | 
| `emos_networks.ipynb` | An implementation of traditional global EMOS using SGD in theano and keras. |
| `fc_and_nn_networks.ipynb` | Experiments with fully connected networks and neural networks. Also contains the experiments with embeddings and auxiliary data. |
| `rnn.ipynb` | Experiments with recurrent neural networks. |
| `categorical_networks.ipynb` | Experiments with categorical/discrete output. |
| `feature_importance.ipynb` | Exploration of feature importance using linear networks and XGBoost. |
| `keras_models.py` | Definition of keras networks. |
| `utils.py` | Helper functions: Data loading, etc. |
| `emos_network_theano.py` | Definition of the EMOS network class in theano. |
| `losses.py` | Definition of custom CRPS loss functions. |

The reference period is 2016, trained with data from 2015. The results are saved in the results directory.