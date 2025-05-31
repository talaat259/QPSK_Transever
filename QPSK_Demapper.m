function [bit1, bit2] = QPSK_Demapper(symbol)
    % Normalize symbol power to match TX normalization
    symbol = symbol * sqrt(2);

    % Decision on bit1 based on real part
    if real(symbol) > 0
        bit1 = 0;
    else
        bit1 = 1;
    end

    % Decision on bit2 based on imaginary part
    if imag(symbol) > 0
        bit2 = 0;
    else
        bit2 = 1;
    end
end
