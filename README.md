# Auxiliary-Function-Based Independent Vector Analysis with Iterative Source Steering (AuxIVA-ISS)

## About
MATLAB script for auxiliary-function-based independent vector analysis with iterative source steering (AuxIVA-ISS) and its application to blind audio source separation.

## Contents
- input [dir]:		includes test audio signals (reverberation time is around 300 ms)
- AuxIVAISS.m:		function of AuxIVA-ISS with pre- and post-processing (STFT, whitening, back projection, and ISTFT)
- backProjection.m:	back projection technique (fixing frequency-wise scales)
- ISTFT.m:			inverse short-time Fourier transform
- main.m:			main script with parameter settings
- STFT.m:			short-time Fourier transform
- whitening.m:		applying principal component analysis for decorrelating observed multichannel signal

## Copyright Note
Coded by Daichi Kitamura
* T. Kim, H. T. Attias, S.-Y. Lee, and T.-W. Lee, "Blind source separation exploiting higher-order frequency dependencies," IEEE Trans. Audio, Speech, and Language Processing, vol. 15, no. 1, pp. 70–79, 2007.
* N. Ono, "Stable and fast update rules for independent vector analysis based on auxiliary function technique", IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA), pp.189-192, 2011.
* S. Robin and N. Ono, "Fast and stable blind source separation with rank-1 updates", IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), pp.236-240, 201j20.

## Python Script
You can find Python script of AuxIVA-IP and AuxIVA-ISS in Pyroomacoustics: https://pyroomacoustics.readthedocs.io/en/pypi-release/pyroomacoustics.bss.auxiva.html#module-pyroomacoustics.bss.auxiva

## See Also
* HP: http://d-kitamura.net