%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%       RP Equation Solver - Mikic model - Verification         %
%       Developed by:   Xinyi Huang                             %
%       Date:           08/21/2017                              %
%       Modified:       08/21/2017                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve the Rayleigh-Plesset equation with Mikic model
% Add global parameters
GLOBALs;

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
R_o0 = 8e-4;

%Tsat = find_sat(water, 'Temperature_sat', pinf*1e-3);
Delta_T = 9.01;
Tinf = 357.23;
Tsat = Tinf - Delta_T;
percentage_water = 0.3;

h_fg = (find_sat(water, 'Enthalpy_gas', Tsat) - find_sat(water, 'Enthalpy_liq', Tsat))*1e3;
rho_gas = find_sat(water,'rho_gas',Tsat);
rho_liq = find_sat(water,'rho_liq',Tsat)*0.3 + find_sat(hexane,'rho_liq',Tsat)*0.7;
rho_water_liq = find_sat(water,'rho_liq',Tsat);
rho_hexane_liq = find_sat(hexane,'rho_liq',Tsat);

thermal_cond = find_sat(water,'Thermal_Conductivity',Tsat)*0.3 + find_sat(hexane,'Thermal_Conductivity',Tsat)*0.7;
capacity = find_sat(water,'Capacity_liq',Tsat)*1e3*0.3 + find_sat(water,'Capacity_liq',Tsat)*1e3*0.7;
thermal_diff = thermal_cond/(rho_liq*capacity);
Ja = rho_liq*capacity*Delta_T/(rho_gas*h_fg);

pi = 3.14159265358;
A = sqrt((2*Delta_T*h_fg*rho_gas) / (3*Tsat*rho_liq));
B = sqrt(12/pi*Ja^2*thermal_diff);

t = [0:1e-6:0.008]';
t_plus = t*A^2/B^2;
R = B^2/A*2/3* (sqrt((t_plus + 1.).^3) - sqrt(t_plus.^3) - 1);

%% Get Perturbation
surface_tension = find_sat(water,'surface_tension',Tsat);
rho = rho_liq;
c = find_sat(water,'sound_speed',Tsat);

R_i = R;
R_o = nthroot(R_i.^3+(R_o0^3*rho_liq - R_i.^3*rho_gas)/rho_water_liq,3);
V_i = zeros(size(R_i,1),1);
V_o = zeros(size(R_o,1),1);
V_i(2:size(R_i,1)) = (diff(R_i')./diff(t'));
V_o(2:size(R_o,1)) = (diff(R_o')./diff(t'));
K = Perturbation(t, R_i);
num = find(K < 5, 1, 'last');
time = t(num);
radius = R_o(num);

mass_bubble = R_i.^3*(rho_gas)*(4*pi/3);
mass_water = mass_bubble + (R_o.^3-R_i.^3)*(rho_liq)*(4*pi/3)*0.3;
mass_hexane = (R_o.^3-R_i.^3)*(rho_liq)*(4*pi/3)*0.7;
mass_all = mass_water + mass_hexane;
figure;plot(t,mass_bubble, t,mass_water, t,mass_hexane, t,mass_all);
xlabel('time(s)');ylabel('mass(kg)');title('Mass variation'); legend('Gas Bubble', 'Water', 'hexane', 'Total mass');
figure;plot(t, K);
figure;plot(t, omega);
