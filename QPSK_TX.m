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
