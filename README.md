# Differentially-private-ANOVA
A repository for an instatiation of differentially private ANOVA testing. To run these files you will need python3 and the python package matplotlib. 

The file `anova.py` contains the main functions used for our differentially private ANOVA testing. The file `csvreader.py` contains the functions necessary to convert the output of our ANOVA test to a plot.

## Running ANOVA

There are three modes of running `anova.py`.  They write to CSV files with very detailed filenames.  Watch out, they overwrite existing files if they are in the directory!

- `python3 anova.py estimate`: Runs the three-group example, using MSE as the estimate of the variance when generating the null distribution.

```
m1 = 0.35
m2 = 0.5
m3 = 0.65
stddev = 0.15
```

- `python3 anova.py realvar`: Runs the three-group example, using `stddev**2` as the real variance when generating the null distribution.

```
m1 = 0.35
m2 = 0.5
m3 = 0.65
stddev = 0.15
```

- `python3 anova.py noisy`: Runs the six-group example, using MSE as the estimate of the variance when generating the null distribution.

```
means_list = [0.5,0.5,0.5,0.6,0.4,0.43]
stddev = 0.2
```

## Plotting values in CSV files

```
python3 csvreader.py ARtest_estimate_1000runs_0.35m1_0.50m2_0.65m3_0.15var_5epsilons_12counts.csv AR_estimate.png 'Differentially private ANOVA, MSE estimate of variance'

python3 csvreader.py ARtest_realvar_1000runs_0.35m1_0.50m2_0.65m3_0.15var_5epsilons_12counts.csv AR_realvar.png 'Differentially private ANOVA, ground truth variance'

python3 csvreader.py ARtest_noisy_1000runs_6means_0.20var_5epsilons_12counts.csv AR_noisy.png 'ANOVA with smaller effect size'
