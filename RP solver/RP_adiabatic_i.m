%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%       RP Equation Solver - RP adiabatic relation - inner r    %
%       Developed by:   Xinyi Huang                             %
%       Date:           08/08/2017                              %
%       Modified:                                               %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function of modified Rayleigh-Plesset equation
% Conditions: infinite boundary radius, with no heat and mass transfer

function dy = RP_adiabatic_i(t,y)
% Initialize the global parameters
GLOBALs;
gamma = Cp/Cv;

% Initialize dy
dy = zeros(2,1);

% Decide the derivatives
dy(1) = y(2);
ratio = 1/nthroot(1+(R_o0/y(1))^3, 3);          % R_o/R_i
%dp = pgas*(R0i/y(1))^(3*gamma) - pinf - 2*surface_tension/y(1)*(1 + ratio) - 4*viscosity*y(2)/y(1)*(1 - ratio^4);
dp = pgas*(R0i/y(1))^(3*gamma) - pinf - 2*surface_tension/y(1)*(1 + ratio) - 4*viscosity*y(2)/y(1)*(1 - ratio^4);
dy(2) = 1/y(1)*(1/(1-ratio)*(1/rho*dp + 0.5*y(2)^2*(1-ratio^4)) - 2*y(2)^2);



