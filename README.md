# Post-processing with Neural Nets

## Data

### Forecasts

**ECMWF forecasts from TIGGE dataset**

https://software.ecmwf.int/wiki/display/TIGGE/Models

- Variables: T2M
- Time period: 2007-2017 (have to check exact dates)
- Forecast initialization time: 00UTC
- Members: 50
- Forecast lead time: 36h and 48h (valid at 12 and 00UTC)
- area: -10E to 30E; 30N to 70N (large part of Europe centered around Germany)
- resolution: 0.5 degrees

### Observations

**DWD SYNOP stations**

- Number of stations: 537
- Variables: T2M (and altitude for each station)
- There is missing observation data!

## Post-processing

### Standard EMOS

Using a rolling training period of the previous 25 days.

| Method | Description |
| ------ | ----------- |
| Global | Use all stations to train one model |
| Local  | To Do       |

The results are in the directory: `standard_postprocessing/preliminary_results`


### Neural Networks

- Minimize CRPS using Stochastic Gradient Descent (or more sophisitcated method in the future)
- In the future: Add more predictor variables
	- Start with more variables at each station, such as wind, humidity (need to think about this in detail.)
	- Add neighborhood information from the ensemble


| Method | Description |
| ------ | ----------- |
| EMOS analog (In production) | A (not yet neural) network mimicking what EMOS does 



