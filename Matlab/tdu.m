% HEADER
clc;
close all
format long

filein = 'msd-050.tdu'; % input file

% universal constants
R_0 = 8.3144598; % [J/(mol*K)] universal gas constant
g_0 = 9.81; % [m/s^2] standard gravity
unitless='[-]';
global debug;
debug = false;
%   === THRUSTER PARAMETERS ===
% IMPORT FROM FILE
fid  = fopen(filein,'r');
if debug;fprintf('reading data from input file (%s)...\n',filein);end
prop_name   = fscanf(fid,'%s',[1,1]); % descriptive header (no quotes, no spaces)
    if debug;fprintf('\tPropellant:\t\t%8s\n',prop_name);end
prop_params = fscanf(fid,'%g',[1 2]); % scan propellant parameters
    k           = prop_params(1,1); % specific heat ratio
    mw          = prop_params(1,2); % molecular weight
    if debug;fprintf('\tk:\t\t%16g\n\tmw:\t\t%16g\n',k,mw);end
total_params= fscanf(fid,'%g',[1 3]); % scan total/stagnation parameters
    T_0         = total_params(1,1); % total temperature
    P_0         = total_params(1,2); % total pressure
    P_b         = total_params(1,3); % back pressure
    if debug;fprintf('\tT_0:\t%16g\n\tP_0:\t%16g\n',T_0,P_0);end
geom_size     = fscanf(fid,'%g',[1 1]); % number of geometry nodes
xcoord = zeros(geom_size,1); radius = zeros(geom_size,1);
    for i=1:geom_size
        geom   = fscanf(fid,'%g',[1 2]);
        xcoord(i)   = geom(1,1); % x coordinate of geometry node
        radius(i)   = geom(1,2); % radius at xcoord
    end
    if debug
        fprintf('\tinlet radius:\t%8f\n\tthroat radius:\t%8f\n\texit radius:\t%8f\n',radius(1),min(radius),radius(end));
        fprintf('\tlength:\t%16f\n\tgeometry nodes:\t%8i\n',xcoord(end),geom_size);
    end
fclose('all'); %close input file
if debug;fprintf('input file closed.\n');end
% MANUAL ENTRY
mw_units = '[kg/mol]';
temperature_units = '[K]';
pressure_units = '[Pa]';
length_units = '[m]';
angle_units = '[deg]';
mdot_units = '[kg/s]';
%   ===------------===   

%           NOZZLE NOMENCLATURE
%   ********************************
%
%              /-
%             /         0 = total parameters
%     ===\---/          1 = chamber
%     (1) (2)  (3)      2 = throat
%     ===/---\          3 = exit
%             \
%              \-
%   ********************************

% % ASSUMPTIONS
% - Uniform flow across nozzle cross sections
% - Isentropic flow (except when normal shock occurs)
% - Initial temperature and pressure are total parameters
% - Inlet flow is unmoving
% 

R = R_0/mw; % gas constant
R_units = '[J/kg K]';
%     A_e = pi*exit_radius^2; % exit area
%     A_t = pi*throat_radius^2; % throat area
A=pi.*radius.^2;
area_units = '[m^2]';

[A_t,A_t_idx]=min(A);
T_star=T_0*(2/(k+1)); %K
P_star=(2/(k+1))^(k/(k-1)); %Pa
rho_star=P_star/(R*T_star); % kg/m^3
mdot=rho_star*sqrt(k*R*T_star)*A_t; %choked

M_idx=linspace(.1,4,length(xcoord));
k1=(2/(k+1));
k2=((k-1)/2);
k3=(.5*(k+1)/(k-1));
% area_ratio=(1./M_idx).*((2/(k+1))*(1+((k-1)/2)*M_idx.^2)).^(.5*(k+1)/(k-1));

if debug;fprintf('Checking for normal shocks...\n');end
M_sub=arearatio2mach_sub(A(end),A_t,k);
M_sup=arearatio2mach_sup(A(end),A_t,k);
pres_ratio_sub=(1+((k-1)/2)*M_sub^2)^(k/(k-1));
pres_ratio_sup=(1+((k-1)/2)*M_sup^2)^(k/(k-1));
Pe1=P_0/pres_ratio_sub; 
Pe2=P_0/pres_ratio_sup;
if P_b >= Pe1
    if debug;fprintf('\tAlways subsonic.\n');end
elseif P_b < Pe1 && P_b > Pe2
    if debug;fprintf('\tNormal shock exists.\n');end
    [ashock, Mminus]=ShockPosition2(A(end),A_t,k,P_0,P_b);
%     xshock=A(find(A==round(ashock*1000)/1000))
%     xshock2=A(find(M==Mminus))
elseif P_b == Pe2
    if debug;fprintf('\tIsentropically supersonic through entire nozzle.\n');end
else
    if debug;fprintf('\tERROR!\n');end
end

if debug;fprintf('Computing flow conditions along x-axis...\n');end
M=zeros(length(xcoord),1);M_sub=M;M_sup=M;temp_ratio=M;T=M;pres_ratio=M;P=M;
M(1)=0;
choked=false;
for x=2:length(xcoord)
        M_sub(x)=arearatio2mach_sub(A(x),A_t,k);
        M_sup(x)=arearatio2mach_sup(A(x),A_t,k);
        if M(x-1)<1
            M(x)=M_sub(x);
        elseif M(x-1)==1
            choked=true;
            if debug;fprintf('\tFlow is choked at A = %f\n',A(x));end
            M(x)=M_sup(x);
        elseif M(x-1)>1
            M(x)=M_sup(x);
        end
        temp_ratio(x)=(1+((k-1)/2)*M(x)^2);
        T(x)=T_0/temp_ratio(x);
        pres_ratio(x)=temp_ratio(x)^(k/(k-1));
        P(x)=P_0/pres_ratio(x); 
        
        temp_ratio_sub(x)=(1+((k-1)/2)*M_sub(x)^2);
        T_sub(x)=T_0/temp_ratio_sub(x);
        pres_ratio_sub(x)=temp_ratio_sub(x)^(k/(k-1));
        P_sub(x)=P_0/pres_ratio_sub(x); 
        temp_ratio_sup(x)=(1+((k-1)/2)*M_sup(x)^2);
        T_sup(x)=T_0/temp_ratio_sup(x);
        pres_ratio_sup(x)=temp_ratio_sup(x)^(k/(k-1));
        P_sup(x)=P_0/pres_ratio_sup(x);
end
thrust = mdot*M(end)*sqrt(k*R*T(end))-(P(end)-P_b)*A(end);
force_units='[N]';
if debug&&not(choked);fprintf('\tFlow is NOT choked.\n');end
if debug;fprintf('done.\n');end
figure
numplots=4;plotcounter=1;
subplot(numplots,1,plotcounter)
plot(xcoord,radius);ylabel('Radius (m)');plotcounter=plotcounter+1;axis([0 xcoord(end) 0 inf]);
%         subplot(numplots,1,plotcounter)
%         semilogy(xcoord,A./A_t,xcoord(mark),A(mark)./A_t,'x');ylabel('A/A_t');plotcounter=plotcounter+1;
subplot(numplots,1,plotcounter)
%         plot(xcoord,M_sub,xcoord,M_sup);ylabel('M');plotcounter=plotcounter+1;axis([0 xcoord(end) 0 inf]);
semilogy(xcoord,M,'k',xcoord,M_sub,xcoord(A_t_idx:end),M_sup(A_t_idx:end));
    legend('M','M<1','M>1');ylabel('M');plotcounter=plotcounter+1;axis([0 xcoord(end) 0 inf]);
subplot(numplots,1,plotcounter)
plot(xcoord,T,'k',xcoord,T_sub,xcoord(A_t_idx:end),T_sup(A_t_idx:end));
    legend('T','T_M_<_1','T_M_>_1');ylabel('T (K)');plotcounter=plotcounter+1;axis([0 xcoord(end) 0 inf]);
%         semilogy(xcoord,T_sub,xcoord,T_sup);ylabel('T');plotcounter=plotcounter+1;axis([0 xcoord(end) 0 inf]);
subplot(numplots,1,plotcounter)
plot(xcoord,P,'k',xcoord,P_sub,xcoord(A_t_idx:end),P_sup(A_t_idx:end));
    legend('P','P_M_<_1','P_M_>_1');ylabel('P (Pa)');plotcounter=plotcounter+1;axis([0 xcoord(end) 0 inf]);
%         semilogy(xcoord,P_sub/10^3,xcoord,P_sup/10^3);ylabel('P');plotcounter=plotcounter+1;axis([0 xcoord(end) 0 inf]);
%         figure
%         semilogy(M(1:find(A==A_t)),A(1:find(A==A_t))./A_t,M(find(A==A_t)+1:end),A(find(A==A_t)+1:end)./A_t,'--');xlabel('M');ylabel('A/A_t');
%         figure
%         semilogy(M_idx,area_ratio,M,A./A_t,'--')
if debug
    figure
        plot(xcoord,radius);ylabel('radius');axis([0 xcoord(end) 0 inf]);
    figure
        semilogy(xcoord,M,'k',xcoord,M_sub,xcoord(A_t_idx:end),M_sup(A_t_idx:end));
        legend('M','M<1','M>1');ylabel('M');axis([0 xcoord(end) 0 inf]);
    figure
        plot(xcoord,T,'k',xcoord,T_sub,xcoord(A_t_idx:end),T_sup(A_t_idx:end));
        legend('T','T_M_<_1','T_M_>_1');ylabel('T (K)');axis([0 xcoord(end) 0 inf]);
    figure
        semilogy(xcoord,P,'k',xcoord,P_sub,xcoord(A_t_idx:end),P_sup(A_t_idx:end),xcoord,ones(length(xcoord))*P_b,'k--');
        legend('P','P_M_<_1','P_M_>_1','P_a_m_b');ylabel('P (Pa)');axis([0 xcoord(end) 0 inf]);
end
%     % format & display outputs
linedivider='------------';
result =  {'Propellant','',prop_name;
           linedivider,'','';
           'Specific heat ratio', k, unitless;
           'Molar mass', mw, mw_units;
           'Specific gas constant',R,R_units;
           linedivider,'','';
           'Total temperature', T_0, temperature_units;
           'Total pressure', P_0, pressure_units;
           linedivider,'','';
           'Length',xcoord(end),length_units;
           'Inlet radius',radius(1),length_units;
           'Throat radius', min(radius), length_units;
           'Exit radius', radius(end), length_units;
           'Inlet area',A(1),area_units;
           'Throat area',A_t,area_units;
           'Exit area',A(end),area_units;
%                'Half-angle',alpha, angle_units;
           linedivider,'','';
           'Throat temperature',T(A_t_idx),temperature_units;
           'Throat pressure',P(A_t_idx),pressure_units;
           'Mass flow rate',mdot,mdot_units;
           linedivider,'','';
           'Exit temperature',T(end),temperature_units;
           'Exit pressure',P(end),pressure_units;
           linedivider,'','';
%                'Exhaust velocity',exit_velocity,velocity_units;
           'Thrust',thrust,force_units;
%                'Specific impulse',specific_impulse,isp_units;
           linedivider,'','';
           'Exit Mach',M(end),unitless;
           'A/At',A(end)/A_t,unitless;
           'T/T0',T(end)/T_0,unitless;
           'P/P0',P(end)/P_0,unitless;
%                'v/at',exit_velocity/throat_velocity,unitless;
           }; 
display(result);
if debug
    fprintf('plotting area ratio...');
    figure
        g=1.4;
        g1=2/(g+1);
        g2=(g-1)/2;
        g3=0.5*(g+1)/(g-1);
        M_span=.0015:.01:7.5;
        Aratio =((1./M_span).*(g1*(1+g2*M_span.^2)).^g3);
        Aratio_sonic =((1./1).*(g1*(1+g2*1.^2)).^g3);
        A_star = A_t/Aratio_sonic; 
        A_theoretical=((1./M_span).*(g1*(1+g2*M_span.^2)).^g3).*A_star;
        A_inlet = ones(length(M_span))*A(1);
        A_exit = ones(length(M_span))*A(end);
        A_throat = ones(length(M_span))*A_t;
        loglog(M_span,A_theoretical,M_span,A_inlet,'--',M_span,A_exit,'--',M_span,A_throat,'--',M_sub,A,'k',M_sup,A,'k')
        legend('theoretical area from F_i_s(M)','inlet area','exit area','throat area','area as function of M found numerically','location','eastoutside')
        grid on
        xlabel('Mach')
        ylabel('Area (m^2)')
    fprintf('done.\n');
end

function Mach = arearatio2mach_sub(A,A_t,k)
% %   solve Mach number from area ratio by Newton-Raphson Method. (assume
% %   subsonic)
% %   https://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/mach.html
%     P = 2/(k+1);
%     Q = 1-P;
%     R = (A/A_t).^2;
%     a = P.^(1/Q);
%     r = (R-1)/(2*a);
%     X = 1/((1+r)+sqrt(r*(r+2)));  % initial guess    
%     diff = 1;  % initalize termination criteria
%     while abs(diff) > .00001
%         F = (P*X+Q).^(1/P)-R*X;
%         dF = (P*X+Q).^((1/P)-1)-R;
%         Xnew = X - F/dF;
%         diff = Xnew - X;
%         X = Xnew;
%     end
%     Mach = sqrt(X);
M=.001:.0005:1;
g=k;
g1=2/(g+1);
g2=(g-1)/2;
g3=0.5*(g+1)/(g-1);
arearatio=1./M.*(g1*(1+g2*M.^2)).^g3;
diff=abs(arearatio-A/A_t);
[res,i]=min(diff);
Mach=M(i);
end
function Mach = arearatio2mach_sup(A,A_t,k)
% %   solve Mach number from area ratio by Newton-Raphson Method. (assume
% %   supersonic)
% %   https://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/mach.html
%     P = 2/(k+1);
%     Q = 1-P;
%     R = (A/A_t).^((2*Q)/P);
%     a = Q.^(1/P);
%     r = (R-1)/(2*a);
%     X = 1/((1+r)+sqrt(r*(r+2)));  % initial guess
%     diff = 1;  % initalize termination criteria
%     while abs(diff) > .00001
%         F = (P*X+Q).^(1/P)-R*X;
%         dF = (P*X+Q).^((1/P)-1)-R;
%         Xnew = X - F/dF;
%         diff = Xnew - X;
%         X = Xnew;
%     end
%     Mach = 1/sqrt(X);
    M=1:.005:12;
    g=k;
    g1=2/(g+1);
    g2=(g-1)/2;
    g3=0.5*(g+1)/(g-1);
    arearatio=1./M.*(g1*(1+g2*M.^2)).^g3;
    diff=abs(arearatio-A/A_t);
    [res,i]=min(diff);
    Mach=M(i);
end
function [shockArea,Mminus]=ShockPosition2(Ae,At,g,P0,Pe)
global debug
% Ae=3;       %exit area
% At=1;       %Throat area
%Me=0.4      %Exit Mach
% P0=1;       %Total P at the inlet
% Pe=0.6;     %Static P at exit
% g=1.4       %isentropic constant
g1=(g-1)/2;
g2=2/(g+1);
g3=0.5*(g+1)/(g-1);

FUN=@(M) 1+g1*M.^2;
FA=@(M) (1./M).*(g2*FUN(M)).^g3;  %isentropic Area/Areat
MSHOCK=@(M) FUN(M)./(g*M.^2-g1);  %shock for M
Pis=@(M) FUN(M).^(g/g-1);
PSHOCK=@(M) 1+2*g/(g+1)*(M.^2-1);

%Arrays of shock areas satisfying upstream and downstream conditions
MM=1:0.01:10;
MMPLUS=MSHOCK(MM);
ASH1=At*FA(MM);

P0SHOCK=(PSHOCK(MM).*Pis(MMPLUS))./Pis(MM);
P02=P0*P0SHOCK;
Me=2/(g-1)*((P02/Pe).^((g-1)/g)-1); Me=sqrt(Me);
ASH2=Ae*FA(Me)./FA(MMPLUS);
if debug
    figure
    plot(MM,ASH1,'--r','LineWidth',3)
    hold on
    plot(MM,ASH2,'--g','LineWidth',3)
    grid on
    xlabel('Mach Number')
    ylabel('A - Shock Area')
    title('Shock Cross Section Area')
end

%Solution as a       min of residual
[RES,i]=min((ASH1-ASH2).^2);
% RES
shockArea=ASH1(i);
Mminus=MM(i);
end
function display(result)
    global debug
    if debug;fprintf('displaying results...\n');end;
    [n,~]=size(result);
    for i = 1:n 
        fprintf('\n%24s\t%15.4g\t%s',result{i,:});
    end
    fprintf('\n')
    if debug;fprintf('done.\n');end;
end