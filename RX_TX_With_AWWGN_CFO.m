%% Parameters
OSR = 4;             % Oversampling rate
fc = 2000;           % Cutoff frequency (Hz)
N_filter = 20;       % Filter order
beta = 0.5;          % Roll-off factor
Rs = 1.98e6;         % Symbol rate
fs = Rs * OSR;       % Sampling frequency

%% Data Generation
N_bits = 1000;        % Number of data bits
bitsPerVector = 2;   % 2 bits per QPSK symbol
n_symbols = N_bits / bitsPerVector;
bits = randi([0, 1], 1, N_bits);

%% Frame Construction with Preamble
final_bit_stream = [];
bitMatrix = reshape(bits, bitsPerVector, n_symbols).';
preamble = Barker_Sequance();
preamble_len = length(preamble);

for i = 1:n_symbols
    bitstream = bitMatrix(i, :);
    stream = [preamble, bitstream];
    final_bit_stream = [final_bit_stream, stream];
end

%% QPSK Modulation
symbols = [];
for i = 1:2:length(final_bit_stream)-1
    bit1 = final_bit_stream(i);
    bit2 = final_bit_stream(i+1);
    symbols(end+1) = QPSK_TX(bit1, bit2);
end

if mod(length(final_bit_stream), 2) == 1
    symbols(end+1) = QPSK_TX(final_bit_stream(end), 0);
end

%% Transmitter Filtering
upfactor = OSR;
filter_coeffs = RRC_FILTER_TX(fs, Rs, N_filter, beta);
upsampled_symbols = upsample(symbols, upfactor);
tx_signal = conv(filter_coeffs, upsampled_symbols);

%% Channel + BER vs SNR Simulation
SNR_range = 0:2:30;
BER_values = zeros(size(SNR_range));

CFO_Hz = 10;
CFO_phase = 2 * pi * CFO_Hz * (0:length(tx_signal)-1) / fs;
cfo_signal = exp(1j * CFO_phase);

for snr_idx = 1:length(SNR_range)
    SNR_dB = SNR_range(snr_idx);

    % Add CFO
    cfo_output = tx_signal .* cfo_signal;

    % Add AWGN
    signal_power = mean(abs(cfo_output).^2);
    noise_power = signal_power / (10^(SNR_dB / 10));
    noise = sqrt(noise_power/2) * (randn(size(cfo_output)) + 1j*randn(size(cfo_output)));
    RX_with_noise = cfo_output + noise;

    %% Receiver: Matched Filtering
    RX_filtered = conv(conj(fliplr(filter_coeffs)), RX_with_noise);
    group_delay = floor(length(filter_coeffs)/2);
    RX_start_sig = RX_filtered(2*group_delay + 1:end - 2*group_delay);

    best_ber = 1;
    best_received_bits = [];

    for offset = 1:upfactor
        downsampled = RX_start_sig(offset:upfactor:end);

        expected_symbols = length(symbols);
        if length(downsampled) > expected_symbols
            downsampled = downsampled(1:expected_symbols);
        end

        bits_at_RX = [];
        for i = 1:length(downsampled)
            S = downsampled(i);
            [bit1, bit2] = QPSK_Demapper(S);
            bits_at_RX = [bits_at_RX, bit1, bit2];
        end

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
            BER = 1;
        end

        if BER < best_ber
            best_ber = BER;
            best_received_bits = no_preamble_bits;
        end
    end

    BER_values(snr_idx) = best_ber;
end

%% Plot BER vs SNR
figure;
semilogy(SNR_range, BER_values, '-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for QPSK with CFO and AWGN');

%% --- Supporting Functions ---
function rrc = RRC_FILTER_TX(fs, Rs, N, beta)
    Ts = 1/Rs;
    T = 1/fs;
    span = (N-1)/2;
    t = (-span:span) * T;
    rrc = zeros(1, N);
    for i = 1:N
        if t(i) == 0
            rrc(i) = (1 + beta*(4/pi - 1));
        elseif abs(t(i)) == Ts/(4*beta)
            rrc(i) = (beta/sqrt(2)) * ((1 + 2/pi)*sin(pi/(4*beta)) + (1 - 2/pi)*cos(pi/(4*beta)));
        else
            numerator = sin(pi*t(i)/Ts*(1 - beta)) + 4*beta*t(i)/Ts * cos(pi*t(i)/Ts*(1 + beta));
            denominator = pi*t(i)/Ts * (1 - (4*beta*t(i)/Ts)^2);
            rrc(i) = numerator / denominator;
        end
    end
    rrc = rrc / sqrt(sum(rrc.^2));
end

function preamble = Barker_Sequance()
    Barker_Seq = [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1];
    preamble = repmat(Barker_Seq, 1, 2);
end

function symbol = QPSK_TX(bit1, bit2)
    if bit1 == 0 && bit2 == 0
        symbol = (1 + 1j)/sqrt(2);
    elseif bit1 == 0 && bit2 == 1
        symbol = (1 - 1j)/sqrt(2);
    elseif bit1 == 1 && bit2 == 1
        symbol = (-1 - 1j)/sqrt(2);
    else
        symbol = (-1 + 1j)/sqrt(2);
    end
end

function [bit1, bit2] = QPSK_Demapper(symbol)
    symbol = symbol * sqrt(2);
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