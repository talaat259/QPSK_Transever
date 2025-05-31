bits = randi([0, 1], 1, 10)
%% 
symbols=[];
for i = 1:2:length(bits)-1
    bit1 = bits(i);
    bit2 = bits(i+1);
    
    % Do something with bit1 and bit2
    fprintf('Pair: %d %d\n', bit1, bit2);
    symbols(end+1)=QPSK_TX(bit1,bit2);
end
symbols;