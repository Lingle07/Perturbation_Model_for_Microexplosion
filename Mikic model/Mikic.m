%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%       RP Equation Solver - Mikic model                        %
%       Developed by:   Xinyi Huang                             %
%       Date:           08/08/2017                              %
%       Modified:       08/17/2017                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%% Solve the Rayleigh-Plesset equation which has viscosity
% Add global parameters
GLOBALs;

% Get water properties: 1-Temperature, 2-Pressure, 3-Density_liquid,
% 12-Viscosity, 14-Surface_Tension, 15-Density_gas
file_water = 'F:\NewlyAdded\UGVR\Rayleigh_Plesset\water_saturate.cgi';
file_hexane = 'F:\NewlyAdded\UGVR\Rayleigh_Plesset\hexane_saturate.cgi';
water = dlmread(file_water,'\t',1,0);
hexane = dlmread(file_hexane,'\t',1,0);

i = 0;
for Temperature = 460:5:600
    i = i+1;
% Decide the ambient temperature(Kelvin) and pressure(Pa)
Tinf = Temperature;
pinf = 1e3;
R_o0 = 1e-4;
evolution_time = 2e-6;

Tsat = find_sat(water, 'Temperature_sat', pinf*1e-3);
Delta_T = Tinf - Tsat;

h_fg = (find_sat(water, 'Enthalpy_gas', Tsat) - find_sat(water, 'Enthalpy_liq', Tsat))*1e3;
rho_gas = find_sat(water,'rho_gas',Tsat);
rho_liq = find_sat(hexane,'rho_liq',Tsat);

thermal_cond = find_sat(water,'Thermal_Conductivity',Tsat);
capacity = find_sat(water,'Capacity_liq',Tsat);
thermal_diff = thermal_cond/(rho_liq*capacity);
Ja = rho_liq*capacity*Delta_T/(rho_gas*h_fg);

pi = 3.14159265358;
A = sqrt((2*Delta_T*h_fg*rho_gas) / (3*Tsat*rho_liq));
B = sqrt(12/pi*Ja^2*thermal_diff);

t = [0:1e-7:1e-4]';
t_plus = t*A^2/B^2;
R = B^2/A*2/3* (sqrt((t_plus + 1.).^3) - sqrt(t_plus.^3) - 1);

%% Get Perturbation
surface_tension = find_sat(water,'surface_tension',Tinf);
rho = rho_liq;
c = find_sat(water,'sound_speed',Tinf);

R_i = R;
R_o = nthroot(R_i.^3+R_o0^3,3);
V_i = zeros(size(R_i,1),1);
V_o = zeros(size(R_o,1),1);
V_i(2:size(R_i,1)) = (diff(R_i')./diff(t'));
V_o(2:size(R_o,1)) = (diff(R_o')./diff(t'));
K = Perturbation(t, R_i);
num = find(K < 5, 1, 'last');
time(i) = t(num);
radius(i) = R_o(num);
end
Temperature = [460:5:600]';
figure;plot(Temperature, time);
figure;plot(Temperature, radius);

% 
% R_i = Y(:,1);
% R_o = nthroot(Y(:,1).^3+R_o0^3,3);
% V_i = Y(:,2);
% V_o = V_i.*(R_i./R_o).^2;
% K = Perturbation(T,Y);
% figure;plot(T, R_i);
% % figure; plot(T,Perturb);
% figure;plot(T,K);
% figure; plot(T,omega);

