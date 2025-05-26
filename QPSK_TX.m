function symbol = QPSK_TX(bit1, bit2)
    if (bit1 == 0 && bit2 == 0)
        symbol = 1 + 1i;
    elseif (bit1 == 0 && bit2 == 1)
        symbol = -1 + 1i;
    elseif (bit1 == 1 && bit2 == 1)
        symbol = -1 - 1i;
    elseif (bit1 == 1 && bit2 == 0)
        symbol = 1 - 1i;
    end

    symbol = symbol / sqrt(2);
end