OSR=4 ;       % Sampling rate (Hz)
fc = 2000;           % Cutoff frequency (Hz)
N = 18;             % Filter order (taps = N+1)
beta = 0.5; 
Rs = 1.98e6;% Roll-off factor
fs = Rs * OSR;

%original number of bits=22400;
N=20;
bits = randi([0, 1], 1, N);
bitsPerVector = 2;
n = length(bits)/bitsPerVector;
final_bit_stream=[];
symbols=[];
bitMatrix = reshape(bits, bitsPerVector, n).';

preamble=Barker_Sequance();
for i=1:n
    bitstream=bitMatrix(i,:);
    disp(size(bitstream));
    stream=[preamble,bitstream];
    final_bit_stream=[final_bit_stream,stream];
end
%%
for i = 1:2:length(final_bit_stream)-1
    bit1 = final_bit_stream(i);
    bit2 = final_bit_stream(i+1);
    
    % Do something with bit1 and bit2
    
    symbols(end+1)=QPSK_TX(bit1,bit2);
end
%%
upfactor=3;

filter_coeffs = RRC_FILTER_TX(fs, Rs, N, beta);
upsampled_symbols=upsample(symbols,upfactor);
final_output=conv(filter_coeffs,upsampled_symbols);
%%

freqz(final_output, 1, 1024, fs);
title('signal after filtration');
%%
% After RRC filtering
plot(real(final_output));
xlabel('Time'); ylabel('Amplitude');
title('Time-Domain Signal After RRC Filtering');
%% channel is here but will be implemented whe  rx produces the correct bitsream
%%
% Matched filtering (Receiver side RRC)
RX_start_sig = conv(fliplr(filter_coeffs),final_output);

% Compensate for group delay introduced by filtering
group_delay = floor(length(filter_coeffs) / 2);
total_delay = 2 * group_delay;
RX_start_sig = RX_start_sig(total_delay + 1:end);

% Downsample to symbol rate
downsampled = RX_start_sig(ceil(upfactor/2):upfactor:end);

% QPSK Demapping
bits_at_RX = [];
for i = 1:length(downsampled)
    S = downsampled(i);
    [bit2, bit1] = QPSK_Demapper(S);  % Make sure ordering matches TX
    bits_at_RX = [bits_at_RX, bit1, bit2];%ana alabt l etnen lw fe moshkla lel mapping
end

% Reshape bits into matrix form: [preamble | data]
preamble_len = 26;
%bitsPerVector = 1;
block_len = preamble_len + bitsPerVector;
expected_rx_bits = n * block_len;

n_bits_matrix = reshape(bits_at_RX(1:expected_rx_bits), block_len, n).';

% Remove preambles and extract data bits
no_preamable_bits = [];
for i = 1:size(n_bits_matrix, 1)
    row = n_bits_matrix(i, :);
    no_preamable_bits = [no_preamable_bits, row(27:end)];
end

% BER Calculation
bit_errors = xor(no_preamable_bits, bits);
disp(['Bit errors: ', num2str(sum(bit_errors))]);
disp(['BER: ', num2str(mean(bit_errors))]);
