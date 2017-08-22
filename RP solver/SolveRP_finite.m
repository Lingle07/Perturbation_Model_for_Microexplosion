%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%       RP Equation Solver                                      %
%       Developed by:   Xinyi Huang                             %
%       Date:           08/08/2017                              %
%       Modified:       08/11/2017                              %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%% Solve the Rayleigh-Plesset equation which has viscosity
% Add global parameters
GLOBALs;

% Get water properties: 1-Temperature, 2-Pressure, 3-Density_liquid,
% 12-Viscosity, 14-Surface_Tension, 15-Density_gas
file_water = 'F:\NewlyAdded\UGVR\Rayleigh_Plesset\water_saturate.cgi';
water = dlmread(file_water,'\t',1,0);

% Decide the ambient temperature(Kelvin) and pressure(Pa)
Tinf = 420;
pinf = 1e5;
R_o0 = 1e-4;
evolution_time = 2e-2;

rho_gas = find_sat(water,'rho_gas',Tinf);
rho_liq = find_sat(water,'rho_liq',Tinf);
pgas = find_sat(water,'pressure',Tinf)*1e3 * rho_liq / rho_gas;
Cp = find_sat(water,'Cp',Tinf);
Cv = find_sat(water,'Cv',Tinf);
c = find_sat(water,'sound_speed',Tinf);

% Decide the density and other things
rho = rho_liq;
surface_tension = find_sat(water,'surface_tension',Tinf);
viscosity = find_sat(water,'viscosity',Tinf)*1e-6;

% Initial conditions
y0 = zeros(2,1);
R0i = 2*surface_tension/(pgas - pinf);
for i = 1:5
    ratio = 1/nthroot((1+(R_o0/R0i)^3),3);
    R0i = 2*surface_tension/(pgas - pinf)*(1+ratio);
end
R0i = R0i*(1+1e-3);
R0o = nthroot(R0i^3+R_o0^3,3);
y0(1) = R0o;
y0(2) = 0.;

% Solve the equation
options = 0;
[T,Y] = ode45(@RP_adiabatic, [0 evolution_time], y0);

R_o = Y(:,1);
R_i = nthroot(Y(:,1).^3-R_o0^3,3);
V_o = Y(:,2);
V_i = V_o.*(R_o./R_i).^2;
K = Perturbation(T,Y);
figure;plot(T, Y(:,1));
% figure; plot(T,Perturb);
figure;plot(T,K);
figure; plot(T,omega);

