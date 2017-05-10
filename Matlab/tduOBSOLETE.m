function varargout = tdu % main function
    clear all;
    clc;
    format long
    %% universal constants
    R_0= 8.3144598; % [J/(mol*K)] universal gas constant
    g_0 = 9.81; % [m/s^2] standard gravity
    unitless='[-]';
    
%%  === EDIT THESE ===
    % gas properties of propellant
    propellant_name = 'Air';
    k = 1.4; % 1.4 for air
    molecular_weight = .0289645; % .0289645 for air
    molar_mass_units = '[kg/mol]';
    
    % chamber conditions
    T_0 = 273; % stagnation temperature
    temperature_units = '[K]';
    P_0 = 101325; % stagnation pressure
    pressure_units = '[Pa]';
    
    % nozzle geometry
    inlet_radius = .0075; % radius at inlet of converging section
    exit_radius = .00708; % radius at exit of diverging section
    throat_radius = .005; % radius at throat
%     exit_radius = .00708; % radius at nozzle exit
%     throat_radius = .005; % radius at nozzle throat
    length_units = '[m]';
    conical_half_angle = 15; % half angle of conical nozzle, 15 degrees is optimal
    angle_units = '[deg]';
%   ===------------===   

%% Math (the fun part)

%           NOZZLE NOMENCLATURE
%   ********************************
%
%              /-
%             /
%     ===\---/          c = chamber
%     (c) (t)  (e)      t = throat
%     ===/---\          e = exit
%             \
%              \-
%   ********************************

    R = R_0/molecular_weight; % gas constant
    R_units = '[J/kg K]';
    A_e = pi*exit_radius^2; % exit area
    A_t = pi*throat_radius^2; % throat area
    area_units = '[m^2]';
    length = (exit_radius-throat_radius)/tan(deg2rad(conical_half_angle)); % [m]
    
    %% Throat conditions
    % We are designing a nozzle that "chokes" the fluid flow at the throat,
    % so we assume that at the throat, Mach number is isentropically
    % brought to the boundary between supersonic and subsonic flow (Mach=1)
    throat_temperature = T_0*(2/(k+1)); 
    throat_pressure = P_0*(2/(k+1))^(k/(k-1)); 
    throat_velocity = sqrt(k*R*throat_temperature);
    % use nozzle areas to find Mach number at the exit via the Newton-Raphson method
    exit_mach = solve_mach(A_e,A_t,k);
%   since
%       mass_flow_rate = density*speed*cross_sectional_area
%   and from the ideal gas law, PV=RT
    throat_density = throat_pressure/(R*throat_temperature); %[kg/m^3]
%   and at the throat Mach = 1 so flow rate (speed of flow) is the fluid's speed of sound at the throat
    throat_flowrate = sqrt(k*R*throat_temperature);
    mass_flowrate = throat_density*throat_flowrate*A_t; 
    mass_flowrate_units = '[kg/s]';

    %% Exit conditions
%   it's much easier to define the flow characteristics at the exit in
%   terms of ratios for now. We'll resolve them into actual useful values
%   later.
%   these are ratios between the throat parameter (numerator) and the exit
%   parameter (denominator). 
%   For example, area_ratio = A_t/A_e
 
%(already know this one)    area_ratio = exit_mach(((k+1)/2)/(1+(k-1)/2*exit_mach^2))^((k+1)/(2*(k-1)));
    temperature_ratio = 1+((k-1)/2)*exit_mach^2;
    pressure_ratio = temperature_ratio^(k/(k-1));
    
    exit_temperature = T_0/temperature_ratio; 
    exit_pressure = P_0/pressure_ratio; 
    
    exit_velocity = sqrt((2*k*R*T_0)/(k-1)*(1-1/(1+(k-1)/2*exit_mach^2)));
    velocity_units = '[m/s]';
    
    %% Thrust and Specific Impulse
    thrust = mass_flowrate*exit_velocity + exit_pressure*A_e;
    force_units = '[N]';
    specific_impulse = thrust/(mass_flowrate*g_0);
    isp_units='[s]';
    
    %% format & display outputs
    linedivider='------------';
    result =  {'Propellant','',propellant_name;
               linedivider,'','';
               'Specific heat ratio', k, unitless;
               'Molar mass', molecular_weight, molar_mass_units;
               'Specific gas constant',R,R_units;
               linedivider,'','';
               'Chamber temperature', T_0, temperature_units;
               'Chamber pressure', P_0, pressure_units;
               'Exit radius', exit_radius, length_units;
               'Throat radius', throat_radius, length_units;
               'Exit area',A_e,'[m^2]';
               'Throat area',A_t,'[m^2]';
               'Half-angle',conical_half_angle, angle_units;
               linedivider,'','';
               'Length',length,length_units;
               'Exit area',A_e,area_units;
               'Throat area',A_t,area_units;
               'Throat temperature',throat_temperature,temperature_units;
               'Throat pressure',throat_pressure,pressure_units;
               'Mass flow rate',mass_flowrate,mass_flowrate_units;
               linedivider,'','';
               'Exit temperature',exit_temperature,temperature_units;
               'Exit pressure',exit_pressure,pressure_units;
               linedivider,'','';
               'Exhaust velocity',exit_velocity,velocity_units;
               'Thrust',thrust,force_units;
               'Specific impulse',specific_impulse,isp_units;
               linedivider,'','';
               'Exit Mach',exit_mach,unitless;
               'A/At',A_e/A_t,unitless;
               'T/Tc',exit_temperature/T_0,unitless;
               'P/Pc',exit_pressure/P_0,unitless;
               'v/at',exit_velocity/throat_velocity,unitless;
               }; 

    display(result);
end

function Mach = solve_mach(A,At,k)
%   solve Mach number from area ratio by Newton-Raphson Method. (assume
%   supersonic)
%   https://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/mach.html
    P = 2/(k+1);
    Q = 1-P;
    R = (A/At).^((2*Q)/P);
    a = Q.^(1/P);
    r = (R-1)/(2*a);
    X = 1/((1+r)+sqrt(r*(r+2)));  % initial guess
    diff = 1;  % initalize termination criteria
    while abs(diff) > .0001
        F = (P*X+Q).^(1/P)-R*X;
        dF = (P*X+Q).^((1/P)-1)-R;
        Xnew = X - F/dF;
        diff = Xnew - X;
        X = Xnew;
    end
    Mach = 1/sqrt(X);
end

function display(result)
    [n,~]=size(result);
    for i = 1:n 
        fprintf('\n%24s\t%15.8f\t%s',result{i,:});
    end
    fprintf('\n')
end
