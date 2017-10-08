# Project Notebook

*This "lab" notebook is for keeping track of our ideas and decisions.*

### October 4

*Stephan*

Yesterday, we (Sebastian, Alexander and I) talked about our results so far. We concluded that the next steps should be:
- Get a standard PP reference that also uses auxiliary variables for a fair comparison with the NN techniques (Sebastian)
- Use XGBoost to estimate the importance of the auxiliary variables (Stephan)
	- Get more auxiliary variables
- Try our convolutional methods

Today I talked to Tom Hammil about which auxiliary variables could be useful for surface temperature PP. He emphasized the importance of the surface energy balance and the variables which are related to that. In the TIGGE dataset there are:
- 2m dewpoint temperature (humidity information)
- Surface winds
- Soil moisture
- Sensible and latent heat flux
- Solar and thermal radiation

There are many more which could be useful. Maybe if the XGBoost method works we can download them all and pick the most important.

### October 8

*Stephan*

Some ideas and concerns:
- We should build a framework for evaluating the model performance that is independent of the post processing scripts. We could define a format for saving the PP predictions such as a CSV file with date, station, mu, sigma. Then one script could read these files, load the observation data and compute the CRPS scores. I believe this would be good to make sure we are not making any mistakes in our scripts! I will start writing something like this soon.
- Training, validation, test split. Currently in my neural network scripts I am only using a training set (2015) and a test set (2016). This is not good machine learning practice. Since I am trying to get the best test score, I will be prone to overfit to my test set. I should therefore not look at my test set until the very end and use a different data set as my validation set. I could eather subdivide 2015 into a training and validation set (but these sets are then probably not independent) or use a different year. For example 2015 for training and 2014 for validation.