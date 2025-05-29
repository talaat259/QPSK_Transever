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
