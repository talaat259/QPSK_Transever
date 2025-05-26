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
