% Example: Call RRCFILTER_TX from another file
 OSR=4 ;       % Sampling rate (Hz)
fc = 2000;           % Cutoff frequency (Hz)
N = 17;             % Filter order (taps = N+1)
beta = 0.5; 
Rs = 1.98e6;% Roll-off factor
fs = Rs * OSR;
% Call the function to get coefficients
filter_coeffs = RRC_FILTER_TX(fs, Rs, N, beta);

% Plot the frequency response
freqz(filter_coeffs, 1, 1024, fs);
title('RRC-Windowed FIR Lowpass Filter');