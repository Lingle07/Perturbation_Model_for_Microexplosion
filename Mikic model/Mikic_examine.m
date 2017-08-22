%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%       RP Equation Solver - Mikic model - Verification         %
%       Developed by:   Xinyi Huang                             %
%       Date:           08/20/2017                              %
%       Modified:       08/21/2017                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve the Rayleigh-Plesset equation which has viscosity
% Get water properties: 1-Temperature, 2-Pressure, 3-Density_liquid,
% 12-Viscosity, 14-Surface_Tension, 15-Density_gas
file_water = 'F:\NewlyAdded\UGVR\Rayleigh_Plesset\water_saturate_1.cgi';
file_hexane = 'F:\NewlyAdded\UGVR\Rayleigh_Plesset\hexane_saturate.cgi';
water = dlmread(file_water,'\t',1,0);
hexane = dlmread(file_hexane,'\t',1,0);

% Decide the ambient temperature(Kelvin) and pressure(Pa)
% B2, B5, B6
% P         1262,   12589.1, 38660.3
% Delta_T   15.74,  10.67,  9.01
% Ja        2690.0, 201.1,  57.5
% Tinf      299.3,  334.24,  357.23
pinf = 38.66e3;
evolution_time = 0.008;

%Tsat = find_sat(water, 'Temperature_sat', pinf*1e-3);
Delta_T = 9.01;
Tinf = 357.23;
Tsat = Tinf - Delta_T;

h_fg = (find_sat(water, 'Enthalpy_gas', Tsat) - find_sat(water, 'Enthalpy_liq', Tsat))*1e3;
rho_gas = find_sat(water,'rho_gas',Tsat);
rho_liq = find_sat(water,'rho_liq',Tsat);

thermal_cond = find_sat(water,'Thermal_Conductivity',Tsat);
capacity = find_sat(water,'Capacity_liq',Tsat)*1e3;
thermal_diff = thermal_cond/(rho_liq*capacity);
Ja = rho_liq*capacity*Delta_T/(rho_gas*h_fg);

pi = 3.14159265358;
A = sqrt((2*Delta_T*h_fg*rho_gas) / (3*Tsat*rho_liq));
B = sqrt(12/pi*Ja^2*thermal_diff);

t = [0:1e-6:0.008]';
t_plus = t*A^2/B^2;
R = B^2/A*2/3* (sqrt((t_plus + 1.).^3) - sqrt(t_plus.^3) - 1);

Examination(:,3) = R;

figure; plot(t, R);
