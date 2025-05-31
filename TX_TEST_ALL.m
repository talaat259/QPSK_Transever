%% Parameters
OSR = 4;             % Oversampling rate
fc = 2000;           % Cutoff frequency (Hz)
N_filter = 20;       % Filter order (renamed to avoid confusion)
beta = 0.5;          % Roll-off factor
Rs = 1.98e6;         % Symbol rate
fs = Rs * OSR;       % Sampling frequency

%% Data Generation
N_bits = 20;         % Number of data bits (renamed for clarity)
bits = randi([0, 1], 1, N_bits);
bitsPerVector = 2;   % 2 bits per QPSK symbol
n_symbols = length(bits) / bitsPerVector;

fprintf('Generated %d bits for %d symbols\n', N_bits, n_symbols);
fprintf('Original bits: %s\n', mat2str(bits));

%% Frame Construction with Preamble
final_bit_stream = [];
bitMatrix = reshape(bits, bitsPerVector, n_symbols).';
preamble = Barker_Sequance();
preamble_len = length(preamble);

fprintf('Preamble length: %d bits\n', preamble_len);

% Build frames: [preamble + data_bits] for each symbol group
for i = 1:n_symbols
    bitstream = bitMatrix(i, :);
    stream = [preamble, bitstream];
    final_bit_stream = [final_bit_stream, stream];
end

fprintf('Total bits in stream: %d\n', length(final_bit_stream));

%% QPSK Modulation
symbols = [];
for i = 1:2:length(final_bit_stream)-1
    bit1 = final_bit_stream(i);
    bit2 = final_bit_stream(i+1);
    symbols(end+1) = QPSK_TX(bit1, bit2);
end

% Handle odd number of bits
if mod(length(final_bit_stream), 2) == 1
    symbols(end+1) = QPSK_TX(final_bit_stream(end), 0);
    fprintf('Warning: Padded with zero for odd bit count\n');
end

fprintf('Number of modulated symbols: %d\n', length(symbols));

%% Transmitter Filtering
upfactor = 3;
filter_coeffs = RRC_FILTER_TX(fs, Rs, N_filter, beta);
upsampled_symbols = upsample(symbols, upfactor);
final_output = conv(filter_coeffs, upsampled_symbols);

fprintf('Filter length: %d\n', length(filter_coeffs));
fprintf('Transmitted signal length: %d\n', length(final_output));

%% Visualization
figure;
freqz(final_output, 1, 1024, fs);
title('Signal After Filtration');

figure;
subplot(2,1,1);
plot(real(final_output));
xlabel('Sample'); ylabel('Amplitude');
title('Real Part - Time-Domain Signal After RRC Filtering');

subplot(2,1,2);
plot(imag(final_output));
xlabel('Sample'); ylabel('Amplitude');
title('Imaginary Part - Time-Domain Signal After RRC Filtering');

%%


%% Receiver Processing with timing offset search

% Matched filtering (Receiver side RRC)
RX_filtered = conv(conj(fliplr(filter_coeffs)), final_output);

fprintf('Received signal length after matched filter: %d\n', length(RX_filtered));

% Compensate for group delay (both TX and RX filter delays)
group_delay = floor(length(filter_coeffs) / 2);
total_delay = 2 * group_delay;  % TX + RX filter delays

if length(RX_filtered) > total_delay
    RX_start_sig = RX_filtered(total_delay + 1:end - total_delay);
else
    RX_start_sig = RX_filtered;
    fprintf('Warning: Signal too short for proper delay compensation\n');
end

% Try all possible downsample offsets to find best BER
best_ber = 1;
best_offset = 1;
best_received_bits = [];

for offset = 1:upfactor
    downsampled = RX_start_sig(offset:upfactor:end);

    % Limit to expected symbols
    expected_symbols = length(symbols);
    if length(downsampled) > expected_symbols
        downsampled = downsampled(1:expected_symbols);
    elseif length(downsampled) < expected_symbols
        fprintf('Warning: fewer symbols than expected at offset %d\n', offset);
    end

    bits_at_RX = [];
    for i = 1:length(downsampled)
        S = downsampled(i);
        [bit1, bit2] = QPSK_Demapper(S);  % Corrected bit order here!
        bits_at_RX = [bits_at_RX, bit1, bit2];
    end

    % Frame processing
    block_len = preamble_len + bitsPerVector;
    expected_rx_bits = n_symbols * block_len;

    if length(bits_at_RX) > expected_rx_bits
        bits_at_RX = bits_at_RX(1:expected_rx_bits);
    elseif length(bits_at_RX) < expected_rx_bits
        bits_at_RX = [bits_at_RX, zeros(1, expected_rx_bits - length(bits_at_RX))];
    end

    n_bits_matrix = reshape(bits_at_RX, block_len, n_symbols).';

    no_preamble_bits = [];
    for i = 1:size(n_bits_matrix, 1)
        row = n_bits_matrix(i, :);
        data_bits = row(preamble_len + 1:end);
        no_preamble_bits = [no_preamble_bits, data_bits];
    end

    if length(no_preamble_bits) == length(bits)
        bit_errors = xor(no_preamble_bits, bits);
        BER = mean(bit_errors);
    else
        BER = 1; % penalize invalid reshapes
    end

    if BER < best_ber
        best_ber = BER;
        best_offset = offset;
        best_received_bits = no_preamble_bits;
    end
end

fprintf('Best downsample offset: %d, BER: %.4f\n', best_offset, best_ber);

% Show final results for best offset
fprintf('\n--- RESULTS (best offset) ---\n');
fprintf('Transmitted bits: %s\n', mat2str(bits));
fprintf('Received bits:    %s\n', mat2str(best_received_bits));
bit_errors = xor(best_received_bits, bits);
num_errors = sum(bit_errors);
fprintf('Bit errors: %d out of %d bits\n', num_errors, length(bits));
fprintf('BER: %.4f (%.2f%%)\n', best_ber, best_ber*100);
if num_errors > 0
    fprintf('Error positions: %s\n', mat2str(find(bit_errors)));
else
    fprintf('Perfect transmission - no errors!\n');
end

%% Supporting functions below

function rrc = RRC_FILTER_TX(fs, Rs, N, beta)
    % fs: sampling rate
    % Rs: symbol rate
    % N: number of taps (must be odd)
    % beta: roll-off factor

    Ts = 1/Rs;                 % Symbol duration
    T = 1/fs;                  % Sampling interval
    span = (N-1)/2;            % Span in samples on either side
    t = (-span:span) * T;      % Time vector

    rrc = zeros(1, N);

    for i = 1:N
        if t(i) == 0
            rrc(i) = (1 + beta*(4/pi - 1));
        elseif abs(t(i)) == Ts/(4*beta)
            rrc(i) = (beta/sqrt(2)) * ((1 + 2/pi)*sin(pi/(4*beta)) + ...
                                       (1 - 2/pi)*cos(pi/(4*beta)));
        else
            numerator = sin(pi*t(i)/Ts*(1 - beta)) + ...
                        4*beta*t(i)/Ts * cos(pi*t(i)/Ts*(1 + beta));
            denominator = pi*t(i)/Ts * (1 - (4*beta*t(i)/Ts)^2);
            rrc(i) = numerator / denominator;
        end
    end

    % Normalize energy
    rrc = rrc / sqrt(sum(rrc.^2));
end

function preamble = Barker_Sequance()
    Barker_Seq = [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1];
    preamble = repmat(Barker_Seq, 1, 2);
end

function symbol = QPSK_TX(bit1, bit2)
    % Map bits to QPSK constellation points (Gray coding)
    % 00 -> (1 + 1j)/sqrt(2)
    % 01 -> (1 - 1j)/sqrt(2)
    % 11 -> (-1 - 1j)/sqrt(2)
    % 10 -> (-1 + 1j)/sqrt(2)

    if bit1 == 0 && bit2 == 0
        symbol = (1 + 1j)/sqrt(2);
    elseif bit1 == 0 && bit2 == 1
        symbol = (1 - 1j)/sqrt(2);
    elseif bit1 == 1 && bit2 == 1
        symbol = (-1 - 1j)/sqrt(2);
    else % bit1==1 && bit2==0
        symbol = (-1 + 1j)/sqrt(2);
    end
end

function [bit1, bit2] = QPSK_Demapper(symbol)
    symbol = symbol * sqrt(2);

    % Use real and imaginary parts for decision
    if real(symbol) > 0
        bit1 = 0;
    else
        bit1 = 1;
    end

    if imag(symbol) > 0
        bit2 = 0;
    else
        bit2 = 1;
    end
end
