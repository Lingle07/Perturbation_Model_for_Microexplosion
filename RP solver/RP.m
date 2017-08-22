%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%       RP Equation Solver - Initial RP function                %
%       Developed by:   Xinyi Huang                             %
%       Date:           08/08/2017                              %
%       Modified:                                               %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function of modified Rayleigh-Plesset equation
% Conditions: infinite boundary radius, with no heat and mass transfer

function dy = RP(t,y)
% Initialize the global parameters
GLOBALs;

% Initialize dy
dy = zeros(2,1);

% Decide the derivatives
dy(1) = y(2);
dy(2) = 1/y(1)*(1/rho*(pgas*(R0/y(1))^3 - pinf - 2*surface_tension/y(1) - 4*viscosity*y(2)/y(1)) - 1.5*y(2)^2);

%dy(1) = 3*y(2);
%dy(2) = 2*nthroot(y(1),3);

