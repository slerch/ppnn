## Overview of results 

Models with values in parentheses are not included in the paper.

### Training 2015 only

| Model        | variant         | CRPS  | computation time (minutes) |
| :--- |:-------------| :-----:| :-----:|
| EMOS  | global          | 1.01 | 0.03 |
|       | (global, window)| (1.00) | (3) |
|       | local           | 0.90 | 0.1 |
|       | (local, window) | (0.90) | (10) |
| EMOS + boosting | local | 0.85 | 14 |
| QRF   | local           | 0.95 |  8 |
| ---   | ---   | --- |  --- |
| Linear network | fc     | 1.01 | 0.2 |
|       | fc_aux          | 0.92 | 0.7 |
|       | fc_emb          | 0.91 | 0.8 |
|       | fc_aux_emb      | 0.88 | 0.8 |
| Neural network | nn_aux_emb | 0.82 | 9 |

### Training 2007-2015

| Model        | variant         | CRPS  | computation time (minutes) |
| :--- |:-------------| :-----:| :-----:|
| EMOS  | global | 1.00 | 0.3 |
|       | (global, window) | (1.00) | (12) |
|       | local | 0.90 | 1 |
|       | (local, window) | (0.88) | (45) |
| EMOS + boosting | local  | 0.80 |  48 |
| QRF   | local   |  0.81 | 430  |
| ---   | ---   | --- |  --- |
| Linear network | fc     | 1.01 | 1 |
|       | fc_aux          | 0.91 | 2 |
|       | fc_emb          | 0.91 | 3 |
|       | fc_aux_emb      | 0.87 | 3 |
| Neural network | nn_aux | 0.87 | 25 |
|  | nn_aux_emb | 0.78 | 16 |
