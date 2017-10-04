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