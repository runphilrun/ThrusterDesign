function varargout = tdu % main function
    clear all;
    clc;
    format long
    % universal constants
    R_0 = 8.3144598; % [J/(mol*K)] universal gas constant
    g_0 = 9.81; % [m/s^2] standard gravity
    unitless='[-]';
    global debug;
    debug = true;
%   === THRUSTER PARAMETERS ===

% % MANUAL ENTRY
% %     % gas properties of propellant
% %     prop_name = 'Air';
% %     k = 1.4; % 1.4 for air
% %     mw = .0289645; % .0289645 for air
% %     mw_units = '[kg/mol]';
% %     
% %     % chamber conditions
% %     T_0 = 273; % stagnation temperature
% %     temperature_units = '[K]';
% %     P_0 = 101325; % stagnation pressure
% %     pressure_units = '[Pa]';
% %     
% %     % nozzle geometry
% %     inlet_radius = .0075; % radius at inlet of converging section
% %     exit_radius = .00708; % radius at exit of diverging section
% %     throat_radius = .005; % radius at throat
% % %     exit_radius = .00708; % radius at nozzle exit
% % %     throat_radius = .005; % radius at nozzle throat
% %     length_units = '[m]';
% %     conical_half_angle = 15; % half angle of conical nozzle, 15 degrees is optimal
% %     angle_units = '[deg]';
%   ===------------===   

%  Math (the fun part)

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

    R = R_0/mw; % gas constant
    R_units = '[J/kg K]';
    A_e = pi*exit_radius^2; % exit area
    A_t = pi*throat_radius^2; % throat area
    area_units = '[m^2]';
    length = (exit_radius-throat_radius)/tan(deg2rad(conical_half_angle)); % [m]
    
    % Throat conditions

    
    % format & display outputs
    linedivider='------------';
    result =  {'Propellant','',prop_name;
               linedivider,'','';
               'Specific heat ratio', k, unitless;
               'Molar mass', mw, mw_units;
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
function loadfromfile(filein)
%  read data from input file
   fid  = fopen(filein,'r');
   %%fprintf('reading data from input file %s...\n',filein);
   header            = fscanf(fid,'%*s',[1 1]) % descriptive header (no quotes, no spaces)
   prop_params       = fscanf(fid,'%s\t%g\t%g\t%s',[1 4]) % scan propellant parameters
    prop_name        = prop_params(1,1) % name of propellant
    k                = prop_params(1,2) % specific heat ratio
    mw               = prop_params(1,3) % molecular weight
%     mw_units
%    
%    g                 = fscanf(ini,'%g',[1 1]);        % body force acceleration (gravity)
%    nodes=zeros(numnod,4);
%    for i=1:numnod
%     nodes(i,:)       = fscanf(ini,'%g',[1 4]);        % node	xcoord  bcflag  bcvalue
%    end
%    %%fprintf('\tnodes imported.\n');
%    elems=zeros(numel,6);
%    for i=1:numel
%     elems(i,:)        = fscanf(ini,'%g',[1 6]);    % elem    node1   node2   dia     E   dens
%    end
   %%fprintf('\telements imported.\n');
   fclose('all'); %close input file
   %fprintf('\tinput file closed.\n');
end
function display(result)
    [n,~]=size(result);
    for i = 1:n 
        fprintf('\n%24s\t%15.8f\t%s',result{i,:});
    end
    fprintf('\n')
end
[ k mw mw_units] = 