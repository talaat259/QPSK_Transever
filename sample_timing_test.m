best_BER = Inf;
best_offset = 1;

for offset = 1:upfactor
    downsampled = RX_start_sig(offset:upfactor:end);

    bits_at_RX = [];
    for i = 1:length(downsampled)
        S = downsampled(i);
        [bit1, bit2] = QPSK_Demapper(S);
        bits_at_RX = [bits_at_RX, bit1, bit2];
    end

    n_bits_matrix = reshape(bits_at_RX(1:expected_rx_bits), block_len, n).';
    no_preamable_bits = [];
    for i2 = 1:size(n_bits_matrix, 1)
        row = n_bits_matrix(i2, :);
        no_preamable_bits = [no_preamable_bits, row(27:end)];
    end

    bit_errors = xor(no_preamable_bits, bits);
    BER = mean(bit_errors);

    if BER < best_BER
        best_BER = BER;
        best_offset = offset;
    end
end

fprintf('Best downsample offset: %d, BER: %f\n', best_offset, best_BER);
