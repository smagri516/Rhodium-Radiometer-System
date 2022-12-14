clear all
clc

% Constants for planks equation
T = 1677; %K
c1 = 3.74177E-16;
c2 = 1.4388E-2;
Area_bot = 2.051017582;

% Calculate f1, Fno, Rlens
Rdet = .0001; %m
Rspot = .3; %m
m = -Rdet/Rspot
z_ml_red = -1.98333333;
z_ld_red = m * z_ml_red;

f1=((1/z_ld_red)-(1/z_ml_red))^-1;

Fno = 4;
Dlens = f1/Fno;
Rlens = Dlens/2;




%% 
y = [0.7
0.71
0.72
0.73
0.74
0.75
0.76
0.77
0.78
0.79
0.8
0.81
0.82
0.83
0.84
0.85
0.86
0.87
0.88
0.89
0.9
0.91
0.92
0.93
0.94
0.95
0.96
0.97
0.98
0.99
1
1.2
1.214
1.223
1.232
1.241
1.251
1.261
1.271
1.281
1.291
1.301
1.312
1.322
1.333
1.344
1.356
1.367
1.379
1.39
1.403
1.415
1.427
1.44
1.453
1.466
1.479
1.493
1.507
1.521
1.536
1.55
1.565
1.573
2
2.012
2.025
2.038
2.051
2.064
2.077
2.09
2.104
2.118
2.132
2.146
2.16
2.174
2.189
2.204
2.219
2.234
2.25
2.266
2.282
2.298
2.314
2.331
2.348
2.365
2.382
2.4
2.418
2.436
2.454
2.473
2.492
2.512
3.146
3.176
3.208
3.273
3.306
3.34
3.375
3.411
3.447
3.484
3.522
3.56
3.6
3.64
3.682
3.724
3.767
3.812
3.857
3.904
3.951
4
4.05
4.101
4.154
4.378
4.438
4.5
4.563
4.628
4.696
4.765
4.836
4.909
4.984
5.062
5.143
5.226
5.311
5.4
5.491];


%% 

y = y';
yopti_append = linspace(.38,.69, 500);
yliwr_append = linspace(5.5,15,1000);
y = [yopti_append y yliwr_append];
y = y.*1E-6; % Convert wavelength to m

%% 
 
e = [0.313
0.31
0.309
0.307
0.306
0.304
0.301
0.303
0.301
0.299
0.298
0.297
0.295
0.292
0.292
0.29
0.29
0.288
0.287
0.285
0.284
0.282
0.282
0.281
0.277
0.277
0.275
0.276
0.274
0.275
0.276
0.262
0.263
0.258
0.258
0.256
0.254
0.252
0.25
0.249
0.247
0.245
0.243
0.24
0.237
0.235
0.236
0.237
0.233
0.233
0.236
0.235
0.235
0.234
0.234
0.234
0.234
0.234
0.233
0.233
0.233
0.233
0.233
0.233
0.243
0.242
0.242
0.243
0.243
0.243
0.243
0.243
0.241
0.241
0.24
0.239
0.239
0.238
0.237
0.236
0.235
0.234
0.234
0.233
0.232
0.231
0.231
0.23
0.23
0.229
0.229
0.229
0.229
0.229
0.229
0.229
0.229
0.229
0.191
0.189
0.185
0.192
0.187
0.182
0.18
0.183
0.184
0.185
0.184
0.183
0.183
0.184
0.181
0.184
0.179
0.179
0.177
0.178
0.176
0.176
0.176
0.176
0.173
0.179
0.171
0.171
0.171
0.171
0.171
.168
0.167
0.166
0.166
0.164
0.164
0.164
0.164
0.164
0.164
];



%% 

e = e';
eopti_append = linspace(.32,.32,500); % Add in epsilon for .38-.7 um
eliwr_append = linspace(.16, .16, 1000);
e = [eopti_append e eliwr_append];

Ly = (c1 ./ (pi .* y.^5)  ) .* (1./ (exp(c2./(y.*T))-1)); % black body
Ly_e = Ly .* e; %gray body

% % Emissivity vs wavelength
% plot(y, e);
% title("Emissivity vs. Wavelength");
% xlabel('\lambda (m)');
% ylabel('$\varepsilon$', 'Interpreter','latex');

% % Ly vs wavelength
% plot(y, Ly);
% title("L_\lambda vs. Wavelength");
% xlabel('\lambda (m)');
% ylabel('L_\lambda (W/m^2*sr*m)');
% legend("T = 1677");

% Ly * e vs wavelength
plot(y, Ly_e);
title("L_\lambda * \epsilon vs. Wavelength");
xlabel('\lambda (m)');
ylabel('L_\lambda*\epsilon_\lambda (W/m^2*sr*m)');

% Maximum Ly part ii
Ly_max = max(Ly)
% Total radiance part iii
L_total = trapz(y, Ly_e)

% Total power part iv
Rs = .3;
Phi = L_total * pi^2 *Rs^2

yopti = y(:,1:509); % .38um --> .78um slicing for all opti wavelengths in array
eopti = e(:,1:509); % slicing for all e at opti lambda
yswir = y(:,551:598); % 1.4um --> 3um slice
eswir = e(:,551:598);
ymwir = y(:,599:902); % 3um-8um
emwir = e(:,599:902);
ylwir = y(:,903:end); %8um --> 15um
elwir = e(:,903:end);

%% 
%part f

coeff = (0.088677 * pi^2 /0.0003305^2) * Rlens^2 * Rdet^2;

Ly_opti =  (c1 ./ (pi .* yopti.^5)  ) .* (1./ (exp(c2./(yopti.*T))-1)); % black body opti spectral radiance
Ly_opti_e = Ly_opti .* eopti; % Real spectral radiance in opti spectrum
L_opti = trapz(yopti, Ly_opti_e) % total radiance of opti band
Phi_det_opti = coeff * L_opti

Ly_swir =  (c1 ./ (pi .* yswir.^5)  ) .* (1./ (exp(c2./(yswir.*T))-1)); % black body opti spectral radiance
Ly_swir_e = Ly_swir .* eswir; % Real spectral radiance in opti spectrum
L_swir = trapz(yswir, Ly_swir_e) % total radiance of opti band
Phi_det_swir = coeff * L_swir

Ly_mwir =  (c1 ./ (pi .* ymwir.^5)  ) .* (1./ (exp(c2./(ymwir.*T))-1)); % black body opti spectral radiance
Ly_mwir_e = Ly_mwir .* emwir; % Real spectral radiance in opti spectrum
L_mwir = trapz(ymwir, Ly_mwir_e) % total radiance of opti band
Phi_det_mwir = coeff * L_mwir

Ly_lwir =  (c1 ./ (pi .* ylwir.^5)  ) .* (1./ (exp(c2./(ylwir.*T))-1)); % black body opti spectral radiance
Ly_lwir_e = Ly_lwir .* elwir; % Real spectral radiance in opti spectrum
L_lwir = trapz(ylwir, Ly_lwir_e) % total radiance of opti band
Phi_det_lwir = coeff * L_lwir
%% 
spectral_lines = [385.652
395.886
413.527
421.114
437.48
];

% Find power emitted in 385_band
y_385_band = y(:,6:14); % band 
Ly_385_band = (c1 ./ (pi .* y_385_band.^5)  ) .* (1./ (exp(c2./(y_385_band.*T))-1)); % black body 385 band spectral radiance
e_385_band = e(:,6:14);
Ly_e_385_band = Ly_385_band .* e_385_band; % accounting emissivity values
L_385_band = trapz(y_385_band, Ly_e_385_band);
Phi_total_bot_385_band = L_385_band * pi * Area_bot; % Power from 385_band from whole bot
Phi_total_det_385_band = L_385_band * coeff

% Find power emitted in 395_band
y_395_band = y(:,22:30);
Ly_395_band = (c1 ./ (pi .* y_395_band.^5)  ) .* (1./ (exp(c2./(y_395_band.*T))-1)); % black body 395 band spectral radiance
e_395_band = e(:,22:30);
Ly_e_395_band = Ly_395_band .* e_395_band; % accounting for emissivity
L_395_band = trapz(y_395_band, Ly_e_395_band);
Phi_total_bot_395_band = L_395_band * pi * Area_bot;
Phi_total_det_395_band = L_395_band * coeff

% Find power emitted in 413_band
y_413_band = y(:,50:58);
Ly_413_band = (c1 ./ (pi .* y_413_band.^5)  ) .* (1./ (exp(c2./(y_413_band.*T))-1)); % black body 413 band spectral radiance
e_413_band = e(:,50:58);
Ly_e_413_band = Ly_413_band .* e_413_band;
L_413_band = trapz(y_413_band, Ly_e_413_band);
Phi_total_bot_413_band = L_413_band * pi * Area_bot;
Phi_total_det_413_band = L_413_band * coeff

% Find power emitted in 421_band
y_421_band = y(:,63:71);
Ly_421_band = (c1 ./ (pi .* y_421_band.^5)  ) .* (1./ (exp(c2./(y_421_band.*T))-1)); % black body 421 band spectral radiance
e_421_band = e(:,63:71);
Ly_e_421_band = Ly_421_band .* e_421_band;
L_421_band = trapz(y_421_band, Ly_e_421_band);
Phi_total_bot_421_band = L_421_band * pi * Area_bot;
Phi_total_det_421_band = L_421_band * coeff

% Find power emitted in 437_band
y_437_band = y(:,89:97);
Ly_437_band = (c1 ./ (pi .* y_437_band.^5)  ) .* (1./ (exp(c2./(y_437_band.*T))-1)); % black body 437 band spectral radiance
e_437_band = e(:,89:97);
Ly_e_437_band = Ly_437_band .* e_437_band;
L_437_band = trapz(y_437_band, Ly_e_437_band);
Phi_total_bot_437_band = L_437_band * pi * Area_bot;
Phi_total_det_437_band = L_437_band * coeff

Phi_det_emission_lines = [Phi_total_det_385_band Phi_total_det_395_band Phi_total_det_413_band Phi_total_det_421_band Phi_total_det_437_band]
Phi_det_emission_lines_total = Phi_total_det_385_band + Phi_total_det_395_band + Phi_total_det_413_band + Phi_total_det_421_band + Phi_total_det_437_band
%% 

