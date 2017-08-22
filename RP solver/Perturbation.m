%% Perturbation calculation
function K = Perturbation(T,Y)
GLOBALs;

% Get the disturbance growth rate
%   All the parameters;
We_o = rho*V_o.^2 .*R_i/surface_tension;
We_i = rho*V_i.^2 .*R_i/surface_tension;
Ma_i = V_i/c;
Delta = R_o./R_i;
phi_o = rho_gas/rho_liq;
phi_i = phi_o;
%   All the parameters of the equation
p1 = Delta - Delta.^2 - phi_o.*Delta;
p2 = (-1 + Delta.^4 + phi_o).*sqrt(We_o);
p3 = 2*(Delta.^2 + 1./Delta.^2);
p4 = -3.*phi_i.*We_i./Ma_i.^2.*Delta.^2;
p = [p1 p2+p1.*3.*sqrt(We_i) p3+p2.*3.*sqrt(We_i)+p4 p3.*3.*sqrt(We_i)];
%   Solve the equations and get perturbation
Omega = zeros(size(p1));
for i = 2:size(R_i,1)
    Omega(i,1) = max(real(roots(p(i,:))));
end
omega = Omega./sqrt(rho_liq.*R_i.^3/surface_tension);
omega(1,1) = 0.;
Perturb = zeros(size(p1));
Perturb = cumtrapz(T,omega);
K = R_o0.*exp(Perturb)./(R_o-R_i);