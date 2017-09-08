#############################################################
## demos for the computational neuroscience course
## Balazs B Ujfalussy - 2017
## http://pattern.wigner.mta.hu/participants/balazs-ujfalussy/

cat('there are two simulations in this demo:\n')
cat('1: membrane potential \n')
cat('2: Hodgkin-Huxley equations \n')
id.demo <- readline('press [1] or [2] to choose between them, or press anything else to quit \n')


if (id.demo == '1'){
	source('./Demos/ions_demo.R', chdir=T)
}

if (id.demo == '2'){
	source('./Demos/HH_demo.R', chdir=T)
}


