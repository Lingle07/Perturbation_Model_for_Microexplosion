%% Find liquid saturation properties from the form
% 1-Temperature, 2-Pressure, 3-Density_liquid,
% 12-Viscosity, 14-Surface_Tension, 15-Density_gas
function out = find_sat(liquid, property, temperature)

if(strcmp(property,'Temperature_sat'))
    line_number = find(liquid(:,2) <= temperature, 1, 'last' );
    lambda = (temperature - liquid(line_number,2))/(liquid(line_number+1,2) - liquid(line_number,2));
    at_temperature = liquid(line_number,:)*(1-lambda) + liquid(line_number+1,:)*lambda;
    out = at_temperature(1);
else
line_number = find(liquid(:,1) <= temperature, 1, 'last' );
if(temperature < min(liquid(:,1)) || temperature >= max(liquid(:,1)))
    out = 0;
    disp('Have problem with getting liquid property!\n');
else
    out = 0;
    lambda = (temperature - liquid(line_number,1))/(liquid(line_number+1,1) - liquid(line_number,1));
    at_temperature = liquid(line_number,:)*(1-lambda) + liquid(line_number+1,:)*lambda;
    if(strcmp(property,'rho_liq'))
        out = at_temperature(3);
    elseif(strcmp(property,'rho_gas'))
        out = at_temperature(15);
    elseif(strcmp(property,'surface_tension'))
        out = at_temperature(14);
    elseif(strcmp(property,'pressure'))
        out = at_temperature(2);
    elseif(strcmp(property,'viscosity'))
        out = at_temperature(12);
    elseif(strcmp(property,'Cp'))
        out = at_temperature(21);
    elseif(strcmp(property,'Cv'))
        out = at_temperature(20);
    elseif(strcmp(property,'sound_speed'))
        out = at_temperature(22);
    elseif(strcmp(property,'Enthalpy_gas'))
        out = at_temperature(18);
    elseif(strcmp(property,'Enthalpy_liq'))
        out = at_temperature(6);
    elseif(strcmp(property,'Thermal_Conductivity'))
        out = at_temperature(13);
    elseif(strcmp(property,'Capacity_liq'))
        out = at_temperature(9);
    end
end
end