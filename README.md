# idegrendszeriModellezes

Modelling and data analysis exercices for the computational neuroscience course
http://cneuro.rmki.kfki.hu/education/neuromodel

1. The modelling exercices

The modelling excercices are centered around 'demos'. Each demo explores a particular neuronal phenomenon. 

To quickly get to the exercices run the corresponding demo (e.g., HH_demo.Rmd) in Rstudio (https://www.rstudio.com)

or the file 'run_demo.R' to run the demos in R (https://www.r-project.org).


Demos are in the 'Demos' folder. Each demo has three separate R files with the following structure:
* xxx_demo.R: interface for the demo. Sets the parameters, runs the simulations and plots the result. 
* xxx_sim.R: defines the simulations.
* xxx_consts.R: defining constants.

Currently there are two demos implemented:
* Nernst equation and membrane potential 
* Hodgkin-Huxley equations


2. Data analysis exercices

Data are in the Data folder, the analysis scripts are in the Analysis folder, and they begin with the same name. To start exploring the data, run the corresponding script in the analysis folder. We have the following data:

* Single neuron somatic current injection experiment from Judit Makara - the original, 50 000 Hz sampling rate
* Single neuron somatic current injection experiment from Judit Makara - reduced, 5 000 Hz sampling rate


