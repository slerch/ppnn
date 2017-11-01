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

### October 17

*Stephan*

**New variables**

Looking at the XGBoost feature importance notebook and from Tom Hammil's suggestions, I would get the following additional variables:

*Surface*
- 2m dewpoint temperature
- 10m U/V wind (also potentially as a target for PP later on)
- convective inhibition
- Soil temperature
- Soil moisture
- Surface latent heat flux
- Surface sensible heat flux
- Total precipitation
- Surface net thermal radiation
- Surface net solar radiation

The feature importance exploration suggests that for the standard deviation, the means of 500hPa geopotential and 850hPa humidity are important. I am not sure why this might be the case. Maybe they are simply indicators of the large-scale weather situation.

**Next steps**

I would like to try building a network which outputs the probability for temperature bins. This way there would be no need to prescribe a distribution. For temperature, which is very normally distributed, I would expect this to be worse, but it would provide a "one-size-fits-all" method for non-normally distributed variables. To test this we could look at wind speed (if these data are available from the stations). I would then also specifically use all 50 members, not the mean and standard deviation. 

Open questions for this are:
- Which metric to use for the categorical output that is comparable to the CRPS? 
- How to chose the bins? What about the first and last bin?


### November 1

*Stephan*

Just a quick update of where things are and what I would like to try:

- Categorical output networks
    - Some things do not fit together yet. For example, the full 50 member ensemble input should be able to get the same results as the parameterized input, but it rellz doesn't. Why not?
    - To explore this I would like to look at the individual steps of the network: 50 member to mean and std to categorical output. How complex does the NN need to be in order to achieve these things individually. 
    - There is also a high sensitivity to hyper-parameters (learning rate, batch size) which I did not see (but also did not test) in the other experiments
    - If wind speed data is available, I would like to try that as an example of a non-Gaussian variable. For this we would also need a reference experiment. 
- RNNs
    - I abandoned these for now since they didn't really work that well. But I think there might still be potential. At some point I would like to revisit this notebook again and go through everything again. Potentially I need to do some more hyper-parameter tuning.
- XGBoost
    - For the parameterized output (mean and std) XGBoost does not work, but we could try it for categorical output. It has the log loss built in, but maybe we could implement the CRPS or at least the RPS.
- Additional data
    - Here I would like to find a standard linear network setup which I can use to get the feature importance. 
    - Should the categorical XGBoost work, that would be a great way to find out about the importance of the features.
    - Then we can try picking the most important features with some data denial experiments.