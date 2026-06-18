% 2026-06 - H2-CIMS ion chemistry simulation - Isaac Standen

% Equations: T.D. Thornberry et al. 2013 & Howard et al. 1972
% Rate constants:
% Howard et al. 1972, Payzant et al. 1972, Feshenfeld 1971

% --------------------------------------------------------------

% (1)  O2+ + H2O + M  -> O2+.H2O + M
% (2)  O2+ + O2 + M  <-> O4+ + M
% (3)  O4+ + H2O      -> O2+.H2O + O2
% (4a) O2+.H2O + H2O  -> H3O+.OH + O2
% (4b)                -> H3O+ + OH + O2
% (5)  H3O+.OH + H2O  -> H3O+.H2O + OH

% --------------------------------------------------------------

clc; clear all;
% Initial conditions

% Rate constants (molec. implicit)
k1   = 2.8e-28;   % cm^6/s
k2f  = 2.6e-30;   % cm^6/s - eqn (2) fwd
k2r  = 2*3.6e-14; % cm^3/s - eqn (2) rvs (sample in air (N2) -> ~2*He k2r)
k3   = 2.2e-09;   % cm^3/s
k4a  = 1.9e-09;   % cm^3/s
k4b  = 3.0e-10;   % cm^3/s
k5   = 3.0e-09;   % cm^3/s

% Timings
dt   = 1e-7; % time step, seconds
t    = 0.01;    % length of sim, seconds
nt   = floor(t/dt); % number of time steps
time = [1:nt];

% Concentration of neutrals
P  = 133.3;      % 1 Torr in Pa
R  = 8.314;      % m3⋅Pa⋅K−1⋅mol−1
kb = R/6.02e23;  % boltzman constant m3PaK-1/molecule
T  = 298;        % K
m  = P/kb/T/1e6; % conc. of air in molecules/cm3
h2o    = 1e-6 * m; % molec/cm^3
o2neut = 0.2 * m;  % molec/cm^3

f_source   = 1/60 * 760; % 1 SLPM/60 seconds * 760 Torr
f_quads    = 1000*1e-5; % 1000 Lps turbo * 1e-5 Torr in quadrupoles
LtoCC      = 1e-3;
o2p_source = 1e7/f_quads *LtoCC; % O2+ from source tube
% Sim. uses one pulse of O2+ at beginning of simulation (not replensished)

% Initialising ion vectors

o2plus = nan(nt, 1); % O2+
o2plus(1) = o2p_source;

o4 = nan(nt, 1); % O4+
o4(1) = 0;

o2_h2o = nan(nt, 1); % O2+.H2O
o2_h2o(1) = 0;

h3o_oh = nan(nt, 1); % H3O+.OH
h3o_oh(1) = 0;

h3o_h2o = nan(nt, 1); % H3O+.H2O
h3o_h2o(1) = 0;

h3o = nan(nt, 1); % H3O+
h3o(1) = 0;

% Model 1 - Euler BWD - with O4+ reverse reaction

for i=2:1:nt
    % O2+
    o2plus(i) = (o2plus(i-1) + (k2r*o4(i-1)*m)*dt)/(1+((k1*h2o*m)+...
        (k2f*o2neut*m))*dt);
    % O4+
    o4(i) = (o4(i-1)+(k2f*(o2plus(i-1)*o2neut*m))*dt)/...
        (1+((k2r*m)+(k3*h2o))*dt);
    % O2+.H2O
    o2_h2o(i) = ...
        (o2_h2o(i-1)+((k1*o2plus(i-1)*h2o*m)+(k3*o4(i-1)*h2o))*dt)/...
        (1+((k4a*h2o)+(k4b*h2o))*dt);
    % H3O+.OH
    h3o_oh(i) = (h3o_oh(i-1)+(k4a*o2_h2o(i-1)*h2o)*dt)/(1+(k5*h2o)*dt);
    % H3O+.H2O
    h3o_h2o(i) = h3o_h2o(i-1)+(k5*h3o_oh(i-1)*h2o)*dt;
    % H3O+
    h3o(i) = h3o(i-1)+(k4b*o2_h2o(i-1)*h2o)*dt;
end

save("ionchem.mat");

% Graphing - Model 1

%clc; load("ionchem.mat");

% Plots
figure(1)

subplot(3,2,1)
plot(time*dt, o2plus, 'b'); hold on;
plot(time*dt, o4, 'g')
plot(time*dt, o2_h2o, 'r')
xlabel("seconds"); ylabel("molec/cm^3")
title("H2-CIMS Ion Funnel Kinetics (with O_4^+ reverse reaction)")
legend("O_2^+", "O_4^+", "O_2^+.H_2O", "Location", "NorthEast");
%yscale log;
grid on; box on; hold off;

subplot(3,2,3)
plot(time*dt, h3o, "b"); hold on;
plot(time*dt, h3o_oh, "g")
plot(time*dt, h3o_h2o, "r")
xlabel("seconds"); ylabel("molec/cm^3")
legend( "H_3O^+", "H_3O^+.OH", "H_3O^+.H_2O", "Location", "NorthEast")
grid on; box on; hold off;

subplot(3,2,5)
plot(time*dt, (h3o+h3o_oh+h3o_h2o), "b"); hold on;
plot(time*dt, (h3o+0.5*h3o_oh+0.5*h3o_h2o), "c")
xlabel("seconds"); ylabel("molec/cm^3")
title(["H_3O^+"])
legend("Total declustering", "50% declustering", "location", "NorthEast")
grid on; box on; hold off;

% Model 2 - Euler BWD - without O4+ reverse reaction

for i=2:1:nt
    % O2+
    o2plus(i) = o2plus(i-1)/(1+((k1*h2o*m)+(k2f*o2neut*m))*dt);
    % O4+
    o4(i) = (o4(i-1)+(k2f*o2plus(i-1)*o2neut*m)*dt)/(1+(k3*h2o)*dt);
    % O2+.H2O
    o2_h2o(i) = ...
        (o2_h2o(i-1)+((k1*o2plus(i-1)*h2o*m)+(k3*o4(i-1)*h2o))*dt)/...
        (1+((k4a*h2o)+(k4b*h2o))*dt);
    % H3O+.OH
    h3o_oh(i) = (h3o_oh(i-1)+(k4a*o2_h2o(i-1)*h2o)*dt)/(1+(k5*h2o)*dt);
    % H3O+.H2O
    h3o_h2o(i) = h3o_h2o(i-1)+(k5*h3o_oh(i-1)*h2o)*dt;
    % H3O+
    h3o(i) = h3o(i-1)+(k4b*o2_h2o(i-1)*h2o)*dt;
end

% Graphing - Model 2

% Plots
figure(1)

subplot(3,2,2)
plot(time*dt, o2plus, 'b'); hold on;
plot(time*dt, o4, 'g')
plot(time*dt, o2_h2o, 'r')
xlabel("seconds"); ylabel("molec/cm^3");
title("H2-CIMS Ion Funnel Kinetics (without O_4^+ reverse reaction)")
legend("O_2^+", "O_4^+", "O_2^+.H_2O", "Location", "NorthEast");
grid on; box on; hold off;

subplot(3,2,4)
plot(time*dt, h3o, "b"); hold on;
plot(time*dt, h3o_oh, "g")
plot(time*dt, h3o_h2o, "r")
xlabel("seconds"); ylabel("molec/cm^3")
legend( "H_3O^+", "H_3O^+.OH", "H_3O^+.H_2O", "Location", "NorthEast")
grid on; box on; hold off;

subplot(3,2,6)
plot(time*dt, (h3o+h3o_oh+h3o_h2o), "b"); hold on;
plot(time*dt, (h3o+0.5*h3o_oh+0.5*h3o_h2o), "c")
xlabel("seconds"); ylabel("molec/cm^3")
title(["H_3O^+"])
legend("Total declustering", "50% declustering", "location", "NorthEast")
grid on; box on; hold off;

% Misc.

% V_beta     = 0.067e06;       % Beta ptcl voltage, V
% V_collison = 60;             % Collision voltage, V
% L_foil     = 3e-2;           % Foil length, m
% rad_tube   = 2e-3;           % Tube radius, m
% ci         = 3.7e10;         % One Curie, /s
% activity   = 15 * 1e-3 * ci; % 15 mCi
% ion_prod   = activity * V_beta/V_collison / ...
%                                  (pi * L_foil * (rad_tube)^2); % 1/cm^3*s
% o2plus     = (ion_prod * 0.26)/760;    % Cracking ratio from Takebe 1972
