# Evaluation of post-processed predictions

The `evaluate_predictions.py` script reads the predictions of the post-processing methods as well as the raw ensemble and computes the mean CRPS. At the moment the prediction have to be for all of 2016. 

The post-processed predictions have to be saved in a csv file in the format

	date, station_id, mean, std

The usage of the script is

	python evaluate_predictions.py --data_dir /Volumes/STICK/data/ppnn_data/ --eval_files emos_network_train_2015_pred_2016.csv

