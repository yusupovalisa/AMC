## AMC: Adaptive Model Combination
### A MATLAB/Octave implementation of this model combination scheme

Author: [Alisa Yusupova](http://www.lancaster.ac.uk/lums/people/alisa-yusupova)<br/>
<!-- E-mail: a(.)yusupova(at)lancaster(.)ac(.)uk<br/> -->
Date:     2021-05-05


#### About

This repository contains: 

- An implementation of the two main elements of the Adaptive Model Combination (AMC) method: 

	- Dynamic linear model with adaptive forgetting (discounting) factor.

	- The fixed share version of the ConfHedge algorithm
	  V. V'yugin and V. Trunov. 
	  [Online aggregation of unbounded losses using shifting experts with confidence](https://doi.org/10.1007/s10994-018-5751-z). 
	  *Machine Learning* (108), 425-444 (2019). 


- Code to reproduce results on simulated data for the corresponding paper.
	- First execute the file `reproduce.m`. 
	This will create a number of CSV files in the two sub-directories `afSims/` and `mcSims/`
	for simulation experiments on adaptive forgetting and model combination, respectively.
	
	- To create the actual figures use `Rscript figures.R` 
