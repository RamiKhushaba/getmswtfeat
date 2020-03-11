Multisignal Wavelet Transform-Based Feature Etxraction
======================================================
getwtfeat is a feature extraction algorithm for ***any kind of signals***, although this was mainly developed for myoelectric, a.k.a, Electromyogram (EMG), signal feature extraction for prostheses control. The algorithm employs the ***coefficients*** of the wavelet transform to extract few features including:

* 1 energy features,
* 2 Variance and std features
* 1 waveform length fetaure
* 1 entropy feature 

You will need to specify the window size, window increment, and the sampling frequency. 

![Alt text](waveletTransform.png?raw=true "fTDD")

As this is a matlab function (adding a python version soon), then usage is really simple, just call this function by submitting the signals matrix (denoted as variable x) as input

	feat = getwtfeat(x,winsize,wininc,samplingFreq)

## Inputs
	x 			columns of signals (rows are samples and column are the signals).
	winsize 	window size.
	wininc		how much to slid the windows by.
	SF			sampling frequency in Hz.

## Outputs
	feat	extracted wavelet features from all nodes/signals
