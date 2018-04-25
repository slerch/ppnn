## Overview of results 

(...) = not in paper

### Training 2015 only

| Model        | variant         | CRPS  | computation time (minutes) |
| :--- |:-------------| :-----:| :-----:|
| EMOS  | global          | 1.01 | 0.03 |
|       | (global, window)| (1.00) | (3) |
|       | local           | 0.90 | 0.85 |
|       | (local, window) | (0.90) | (10) |
| EMOS + boosting | local | 0.85 | 14 |
| QRF   | local           | 0.95 |  28 |
| ---   | ---   | --- |  --- |
| Linear network | fc     | 1.01 | |
|       | fc_aux          | 0.92 | |
|       | fc_emb          | 0.91 | |
|       | fc_aux_emb      | 0.87 | |
| Neural network | nn_aux_emb | 0.82 | |

### Training 2007-2015

| Model        | variant         | CRPS  | computation time (minutes) |
| :--- |:-------------| :-----:| :-----:|
| EMOS  | global | 1.00 | 0.3 |
|       | (global, window) | (1.00) | (12) |
|       | local | 0.90 | 1 |
|       | (local, window) | (0.88) | (45) |
| EMOS + boosting | local  | 0.82 |  48 |
| QRF   | local   |   |   |