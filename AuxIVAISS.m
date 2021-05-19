function [estSig, cost] = AuxIVAISS(mixSig,nSrc,sampFreq,fftSize,shiftSize,windowType,nIter,refMic,applyWhitening,drawConv)
% Blind source separation with auxiliary-function-based independent vector analysis (IVA) based on iterative source steering (AuxIVA-ISS)
%
% Coded by D. Kitamura (d-kitamura@ieee.org)
%
% # Original IVA paper
% T. Kim, H. T. Attias, S.-Y. Lee and T.-W. Lee, 
% "Blind source separation exploiting higher-order frequency dependencies,"
% IEEE Trans. ASLP, vol. 15, no. 1, pp. 70–79, 2007.
%
% # Original AuxIVA paper
% N. Ono, 
% "Stable and fast update rules for independent vector analysis based on auxiliary function technique," 
% Proc. WASPAA, pp. 189–192, 2011.
%
% # Original AuxIVA-ISS paper
% S. Robin and N. Ono, 
% "Fast and stable blind source separation with rank-1 updates," 
% Proc. ICASSP, pp.236–240, 2020.
%
% [syntax]
%   [estSig,cost] = AuxIVAISS(mixSig,nSrc,sampFreq,fftSize,shiftSize,windowType,nIter,refMic,applyWhitening,drawConv)
%
% [inputs]
%         mixSig: observed mixture (sigLen x nCh)
%           nSrc: number of sources in the mixture (scalar)
%       sampFreq: sampling frequency [Hz] of mixSig (scalar)
%        fftSize: window length [points] in STFT (scalar, default: next higher power of 2 that exceeds 0.256*sampFreq)
%      shiftSize: shift length [points] in STFT (scalar, default: fftSize/2)
%     windowType: window function used in STFT (name of window function, default: 'hamming')
%          nIter: number of iterations in the parameter update in AuxIVA-ISS (scalar, default: 100)
%         refMic: reference microphone for applying back projection (default: 1)
% applyWhitening: apply whitening to the observed multichannel spectrograms or not (true or false, default: true)
%       drawConv: plot cost function values in each iteration or not (true or false, default: false)
%
% [outputs]
%         estSig: estimated signals (sigLen x nCh x nSrc)
%           cost: convergence behavior of cost function in IVA (nIter+1 x 1)
%

% Arguments check and set default values
arguments
    mixSig (:,:) double
    nSrc (1,1) double {mustBeInteger(nSrc)}
    sampFreq (1,1) double
    fftSize (1,1) double {mustBeInteger(fftSize)} = 2^nextpow2(0.256*sampFreq)
    shiftSize (1,1) double {mustBeInteger(shiftSize)} = fftSize/2
    windowType char {mustBeMember(windowType,{'hamming','hann','rectangular','blackman','sine'})} = 'hamming'
    nIter (1,1) double {mustBeInteger(nIter)} = 100
    refMic (1,1) double {mustBeInteger(refMic)} = 1
    applyWhitening (1,1) logical = true
    drawConv (1,1) logical = false
end

% Error check
[sigLen, nCh] = size(mixSig); % sigLen: signal length, nCh: number of channels
if sigLen < nCh; error("The size of mixSig might be wrong.\n"); end
if nCh < nSrc || nSrc < 2; error("The number of channels must be equal to or grater than the number of sources in the mixture.\n"); end
if fftSize < 1; error("The FFT length in STFT (fftSize) must be a positive integer value.\n"); end
if shiftSize < 1; error("The shift length in STFT (shiftSize) must be a positive integer value.\n"); end
if nIter < 1; error("The number of iterations (nIter) must be a positive integer value.\n"); end
if refMic < 1 || refMic > nCh; error("The reference microphone must be an integer between 1 and nCh.\n"); end

% Apply multichannel short-time Fourier transform (STFT)
[mixSpecgram, windowInStft] = STFT(mixSig, fftSize, shiftSize, windowType);

% Apply whitening (decorrelate X so that the correlation matrix becomes an identity matrix) based on principal component analysis
if applyWhitening
    inputMixSpecgram = whitening(mixSpecgram, nSrc); % apply whitening, where dimension is reduced from nCh to nSrc when nSrc < nCh
else
    inputMixSpecgram = mixSpecgram(:,:,1:nSrc); % when nSrc < nCh, only mixSpecgram(:,:,1:nSrc) is input to IVA so that the number of microphones equals to the number of sources (determined condition)
end

% Apply AuxIVA-ISS
[estSpecgram, cost] = local_AuxIVAISS(inputMixSpecgram, nIter, drawConv);

% Apply back projection (fix scale ambiguity using reference microphone channel)
scaleFixedSepSpecgram = backProjection(estSpecgram, mixSpecgram(:,:,refMic)); % scale-fixed estimated signal

% Inverse STFT for each source
estSig = ISTFT(scaleFixedSepSpecgram, shiftSize, windowInStft, sigLen);
end

%% Local function for AuxIVA-ISS
function [outY, cost] = local_AuxIVAISS(inX, nIter, drawConv)
% [inputs]
%            inX: observed multichannel spectrograms (I x J x M)
%          nIter: number of iterations of the parameter updates
%       drawConv: plot cost function values in each iteration or not (true or false)
%
% [outputs]
%           outY: estimated spectrograms of sources (I x J x N)
%           cost: convergence behavior of cost function in ILRMA (nIter+1 x 1)
%
% [scalars]
%              I: number of frequency bins
%              J: number of time frames
%              M: number of channels (microphones)
%              N: number of sources (equals to M in IVA-based blind source separation)
%
% [matrices]
%            inX: observed multichannel spectrograms (I x J x M)
%              X: permuted observed multichannel spectrograms (M x J x I)
%              W: frequency-wise demixing matrix (N x M x I)
%           outY: estimated multisource spectrograms (I x J x N)
%              Y: permuted estimated multisource spectrograms (N x J x I)
%              P: estimated multisource power spectrograms (N x J x I)
%              R: L2 norms of frequency vector for each source and each time of the power spectrograms (N x J)
%              v: variable vector for iterative source steering (N x 1)
%

% Initialization
[I, J, M] = size(inX); % I:frequency bins, J: time frames, M: channels
X = permute(inX, [3,2,1]); % permuted X whose dimensions are M x J x I
N = M; % N: number of sources, which equals to M in IVA
W = zeros(N, M, I); % frequency-wise demixing matrix
for i = 1:I
    W(:,:,i) = eye(N); % initial demixing matrices are set to identity matrices
end
Y = X; % initial permuted estimated multisource spectrograms (N x J x I)
P = abs(Y).^2; % power spectrograms of Y (N x J x I)
R = sqrt(sum(P, 3)); % N x J, L2 norm of y (Laplace distribution), where R is common for all freq bins
v = zeros(N,1); % v vector for iterative source steering
cost = zeros(nIter+1, 1); % cost function values

% Calculate initial cost function value
if drawConv
    cost(1,1) = local_calcCostFunction(P, W, I, J);
end

% Optimize parameters in AuxIVA-ISS
fprintf('Iteration:    ');
for iIter = 1:nIter
    fprintf('\b\b\b\b%4d', iIter);
    
    %%%%% Update parameters %%%%%
    for n = 1:N
        YY = Y .* conj(Y(n,:,:)); % N x J x I, using implicit expansion (NxJxI .* 1xJxI)
        for i = 1:I
            for nn = 1:N % calculate v vector BEGIN
                d = sum((1./(2*R(nn,:))).*real(YY(n,:,i)), 2) / J; % scalar
                if nn ~= n
                    u = sum((1./(2*R(nn,:))).*YY(nn,:,i), 2) / J; % scalar
                    v(nn,1) = u / d;
                else
                    v(nn,1) = 1 - 1/sqrt(d);
                end
            end % calculate v vector END
            Y(:,:,i) = Y(:,:,i) - v.*Y(n,:,i); % update Y, usnig implicit expansion (Nx1 .* 1xJ)
        end
    end
    P = abs(Y).^2; % power spectrograms of Y (N x J x I)
    R = sqrt(sum(P, 3)); % N x J, L2 norm of y (Laplace distribution), where R is common for all freq bins
    
    %%%%% Calculate cost function value %%%%%
    if drawConv
        for i = 1:I
            W(:,:,i) = Y(:,:,i)*X(:,:,i)' / (X(:,:,i)*X(:,:,i)'); % derived by "Y(:,:,i) = W(:,:,i)*X(:,:,i)"
        end
        cost(iIter+1,1) = local_calcCostFunction(P, W, I, J);
    end
end
outY = permute(Y, [3,2,1]); % permute (N x J x I -> I x J x N)

% Draw convergence behavior
if drawConv
    figure; plot((0:nIter), cost);
    set(gca, 'FontName', 'Times', 'FontSize', 16);
    xlabel('Number of iterations', 'FontName', 'Arial', 'FontSize', 16);
    ylabel('Value of cost function', 'FontName', 'Arial', 'FontSize', 16);
end
fprintf(' AuxIVA-ISS done.\n');
end

%% Local function for calculating cost function value in IVA
function [ cost ] = local_calcCostFunction(P, W, I, J)
logDetAbsW = zeros(I,1);
for i = 1:I
    logDetAbsW(i,1) = log(max(abs(det(W(:,:,i))), eps));
end
cost = sum(sqrt(sum(P, 3)), 'all')/J - 2*sum(logDetAbsW, 'all');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%