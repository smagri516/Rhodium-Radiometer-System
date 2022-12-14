
%THIS MUST BE RUN AFTER 'q2.m'!!!!

%The band that we will be working in is 380-450nm given the
% spectral lines of rhodium. 
% We will be using the Hammamatsu

y_working = y(:,1:531); % All y values from 380nm to 1um "what detector sees"

%Responsivity in positive linear region
y_R1 = y_working(:,1:244);
slope1 = (20-8)/(500E-9 - 380E-9);
R1 = slope1.*(y_R1 - 380E-9) + 8; % Responsivities for 380-530nm
%Responsivity in top region 530-700nm
y_R2 = y_working(:,245:500);
R2=23.5.*linspace(1,1,256);
%Responsivity in neg linear region 700-1000nm
y_R3 = y_working(:,501:531);
slope3 = (2-21)/(1000E-9 - 700E-9);
R3 = slope3.*(y_R3 - 700E-9) + 21;


R_working = [R1 R2 R3]; % Total approx Responsivity 380-1000nm
% plot(y_working, R_working)

Ly_working_band = (c1 ./ (pi .* y_working.^5)  ) .* (1./ (exp(c2./(y_working.*T))-1)); % black body working band spectral radiance
e_working_band = e(:,1:531);
Ly_e_working_band = Ly_working_band .* e_working_band;
L_working_band = trapz(y_working, Ly_e_working_band);
Phi_total_det_working_band = L_working_band * coeff
RxL = R_working .* Ly_e(:,1:531);

i_working = trapz(y_working, RxL) .* coeff % Total signal (amps) expected
i_emission_spectrum = (trapz(y_385_band, RxL(:,6:14)) + trapz(y_395_band, RxL(:,22:30)) + ...
    trapz(y_413_band, RxL(:,50:58)) + trapz(y_421_band, RxL(:,63:71)) + trapz(y_437_band, RxL(:,89:97))) .* coeff
plot(y_working, RxL)

