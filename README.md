# Simulator for Finite-Alphabet Equalization
(c) 2020 Oscar Castañeda and Christoph Studer
e-mail: oc66@cornell.edu & studer@cornell.edu

More information about our research can be found at [http://vip.ece.cornell.edu] and [https://sites.google.com/site/durisi].

### Important information

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

Oscar Castañeda, Sven Jacobsson, Giuseppe Durisi, Tom Goldstein, and Christoph Studer, "Finite-Alphabet MMSE Equalization for All-Digital Massive MU-MIMO mmWave Communication," IEEE Journal on Selected Areas in Communications (J-SAC), vol. 38, no. 9, September 2020, pp. 2128-2141

and clearly mention this in your paper.

### How to start a simulation:

Simply run

```sh
fa_equalizer_sim
```

which starts a simulation of a 256 BS antenna, 16 user, 16-QAM massive MIMO system using two of the finite-alphabet equalizers proposed in our paper, FL-MMSE and FAME-FBS, both using the 1-bit alphabet. As a baseline, the simulator also includes the conventional L-MMSE equalizer, whose entries are not quantized. The simulation runs in an i.i.d. Rayleigh-fading channel. The code also includes two other finite-alphabet equalizers proposed in our paper: FAME-EXH and FAME-SDR. However, note that the code provided for these last two equalizers only works for 1-bit alphabets.

The simulator runs with predefined parameters. You can specify your own system and simulation parameters by passing your own "par"-structure (see the simulator for an example). Note that we use default parameters for the considered system configuration; if you want to run the simulation with different parameters (including the number of FAME-FBS iterations), then please refer to the MATLAB code for other parameter settings.

We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator.

### Notes

* To use the 'FAME_SDR_1b' equalizer, you need to first install CVX [http://cvxr.com/cvx/]

* To reproduce the results from our paper, aside from using the provided parameters, you need to use the following number of trials on par.trials:
	* For the B=256, U=16, 16-QAM system, use par.trials = 1e4
	* For the B=16, U=4, 16-QAM system, use par.trials = 1e3. You also need to run the simulation 10 times, each time with a different par.runId, from 0 to 9 (these simulations can be run in parallel), and then average the error-rate results.
	* For the B=8, U=2, QPSK system, use par.trials = 2e3. You also need to run the simulation 10 times, each time with a different par.runId, from 0 to 9 (these simulations can be run in parallel), and then average the error-rate results.

### Version history
* Version 0.1: oc66@cornell.edu - initial version for GitHub release
