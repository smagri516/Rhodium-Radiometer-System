% Johnson Noise
kb = 1.38E-23;
h = 6.626E-34;
T = 300; % K
R = 10^6;
c = 3E8;
freq_working = c./y_working; % for whole 380-1000nm
B = 1000E-9 - 380E-9;
% freq_working = c./y_working(:,1:116); % for 380-450 nm
% B = 450E-9 - 380E-9;
i_dark = 50E-12;

ij_integrand = freq_working.^2 ./ (exp(h.*freq_working ./ (kb.*T))-1);
ij_integral = trapz(freq_working, ij_integrand);
i_j = (abs(4 * h * ij_integral/R))^.5;

% SHOT Noise
i_shot = (2 * 1.6E-19 * (i_working + i_dark) * B)^.5;

i_noise = (i_shot^2 + i_j^2 + i_dark^2)^.5;