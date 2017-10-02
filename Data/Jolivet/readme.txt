## Datasets for the computational neuroscience course
## 
## Jolivet:
## 	

# Jolivet_stim.RData contains a single matrix 'stim' with 32 rows, each being the injected current (in pA) of a single experiment
# Jolivet_resp.RData contains a single matrix 'resp', with 32 rows, each being the voltage recorded at the soma (mV)

## dataset: 32 traces, 8 different current injection, 4 repetitions each
##			each repetition is 34 s long, sampled with 1000 Hz
## Parameters of the input current were
# m.I <- c(	645,	482,	161,	480,	240,	481,	164,	241) # mean, pA
# sd.I <- c(	328,	164,	164,	42,	41,	6,	328,	6) # sd
# rate <- c(	27,	19,	1,	16,	2,	15,	6,	3) # firing rate
# 		1-4	5-8	9-12,	13-16,	17-20,	21-24	25-28	29-32

