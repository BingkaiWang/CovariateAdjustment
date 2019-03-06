# Covariate Adjustment
This repo is used to store the code for paper "Analysis of Covariance (ANCOVA) in Randomized Trials: More Precision and Valid Confidence Intervals, Without Model Assumptions".

The Data_Preprocessing_and_Analysis folder contains code of raw data preprocessing and data analysis for the MCI, METS and TADS trial.

The Simulation folder contain code for simulation and visualization.

To apply different methods to estimate the average treatment effect, one can use the "adjust_estimator" function, with description given in "Data_Preprocessing_and_Analysis/adjust_estimator.R". An example of using this function is in "example.R". For handling missing data by Doubly-robust estimator, one can refer to the example given in "Data_Preprocessing_and_Analysis/Missing_data.R". 

Since the data from the trials is not public, we cannot include it in this Repository. Instead, we provide code to generate sample data sets in "example.R". Running that file will create a 500-participant data set, apply the unadjusted and ANCOVA estimators, and estimate the sample size reduction from adjustment.

The draft of this paper is currently stored [here](https://biostats.bepress.com/jhubiostat/paper292/), where main results and data analysis of 3 clinical trials are given.

The supporting information of this paper is currently stored [here](http://people.csail.mit.edu/mrosenblum/ANCOVA.pdf), where proofs and formulas are given.
