%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample program for blind source separation based on AuxIVA-ISS          %
%                                                                         %
% Coded by D. Kitamura (d-kitamura@ieee.org)                              %
%                                                                         %
% # Original paper (AuxIVA-IP)                                            %
% N. Ono, "Stable and fast update rules for independent vector analysis   %
% based on auxiliary function technique," Proc. WASPAA, pp.189–192, 2011. %
%                                                                         %
% # Original paper                                                        %
% S. Robin and N. Ono, "Fast and stable blind source separation with      %
% rank-1 updates," Proc. ICASSP, pp.236-240, 2020.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

% Set parameters
seed = 1; % pseudo random seed
refMic = 1; % reference microphone for back projection
resampFreq = 16000; % resampling frequency [Hz]
nSrc = 2; % number of sources
fftSize = 4096; % window length in STFT [points]
shiftSize = 2048; % shift length in STFT [points]
windowType = "hamming"; % window function used in STFT
nIter = 30; % number of iterations (define by checking convergence behavior with drawConv=true)
applyWhitening = false; % true or false (true: apply whitening to the observed multichannel spectrograms)
drawConv = true; % true or false (true: plot cost function values in each iteration and show convergence behavior, false: faster and do not plot cost function values)

% Fix random seed
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed))

% Input data and resample
[srcSig(:,:,1), sampFreq] = audioread('./input/drums.wav'); % signal x channel x source (source image)
[srcSig(:,:,2), sampFreq] = audioread('./input/piano.wav'); % signal x channel x source (source image)
srcSigResample(:,:,1) = resample(srcSig(:,:,1), resampFreq, sampFreq, 100); % resampling for reducing computational cost
srcSigResample(:,:,2) = resample(srcSig(:,:,2), resampFreq, sampFreq, 100); % resampling for reducing computational cost

% Mix source images of each channel to produce observed mixture signal
mixSig(:,1) = srcSigResample(:,1,1) + srcSigResample(:,1,2);
mixSig(:,2) = srcSigResample(:,2,1) + srcSigResample(:,2,2);
if abs(max(max(mixSig))) > 1 % check clipping
    error('Cliping detected while mixing.\n');
end

% Reference signals
refSig(:,1) = srcSigResample(:,refMic,1);
refSig(:,2) = srcSigResample(:,refMic,2);

% Blind source separation based on AuxIVA-ISS
[estSig, cost] = AuxIVAISS(mixSig, nSrc, resampFreq, fftSize, shiftSize, windowType, nIter, refMic, applyWhitening, drawConv);

% Output separated signals
outputDir = sprintf('./output');
if ~isfolder( outputDir )
    mkdir( outputDir );
end
audiowrite(sprintf('%s/observedMixture.wav', outputDir), mixSig(:,refMic), resampFreq); % observed signal
audiowrite(sprintf('%s/originalSource1.wav', outputDir), refSig(:,1), resampFreq); % source signal 1
audiowrite(sprintf('%s/originalSource2.wav', outputDir), refSig(:,2), resampFreq); % source signal 2
audiowrite(sprintf('%s/estimatedSignal1.wav', outputDir), estSig(:,1), resampFreq); % estimated signal 1
audiowrite(sprintf('%s/estimatedSignal2.wav', outputDir), estSig(:,2), resampFreq); % estimated signal 2

fprintf('The files are saved in "./output".\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%