## Overview of results of benchmark models

### Training 2015 only

| Model        | variant         | CRPS  | computation time (minutes) |
| :--- |:-------------| :-----:| :-----:|
| EMOS  | global | 1.01 | 0.75 |
|       | global, window | 1.00 | 3 |
|       | local | 0.91 | 0.85 |
|       | local, window | 0.90 | 10 |
| EMOS + boosting | local  | 0.87 | 14 |
| QRF   | local   | 0.95 |  28 |

### Training 2007-2015

| Model        | variant         | CRPS  | computation time (minutes) |
| :--- |:-------------| :-----:| :-----:|
| EMOS  | global | 1.00 | 1 |
|       | global, window | 1.00 | 12 |
|       | local | 0.90 | 1 |
|       | local, window | 0.88 | 45 |
| EMOS + boosting | local  | 0.82 |  48 |
| QRF   | local   | x | x |