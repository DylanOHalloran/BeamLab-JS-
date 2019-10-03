        %%%%%%%%%%%%%%%% Theo&Exp Free response %%%%%%%%%%%%%%
        %%%%%%%%%%%% 3B5 Vibrating Beam Laboratory %%%%%%%%%%%
        %%%%%%% If in trouble, contact the demonstrator %%%%%%
        %%%%%%%%%% Federico Caruso (fcaruso@tcd.ie) %%%%%%%%%%
clc
clear
close all

fs = 1000;         % sample rate
%% Input data from the lab sheet
M = ; % Concentrated mass, kg
m_cl_ds = ; % Mass of clamps for spring and dashpot, kg
m_motor = ; % Mass of the motor (including clamps and rotor), kg
m_rotor = ; % Mass of the rotor, kg
m_beam = ; % Mass of the beam, kg
l_1 = ; % Length l1, m
l_2 = ; % Length l2, m
l_3 = ; % Length l3, m
l_4 = ; % Length l4, m
l_beam = ; % Length of the beam, m
k = ; % Spring constant, N/m
I_beam = ; % Inertia moment of the beam kg*m^2
f_f = ;     % Measured Frequency of rotating eccentric mass (RPS)
xi = ;     % damping ratio measured from free response
r_rotor = ;   %length of rotor
d_a = ; % distance of the rotor pivot from edge

w_rad = ; % Rotating eccentric mass frequency, rad/s

d_rotor = ;  %distance of rotor cg to rotation axis; %please do not edit this command
t_end = ;    % Forced response sampling time, s;

% Measurement factor
measurement_factor = l_1/l_3; % please do not edit this command

%% Calculate Experimental response

% Load Forced experimental data
forced_exp_res = load(Forced_Response_20th);

%% Separate into two vectors
% Displacement
x_for_res = (forced_exp_res(:,2));
% Time
t_exp = forced_exp_res(:,1);
% Transfer experimental displacements to x1 position and convert m to mm
x_for_res = x_for_res*measurement_factor*0.001;

%% Theoretical forced damped response
% Forced response natural frequency
w_n_rad = 

%% Motor parameters
% Inertia moment about G
I_G = ; %units
% Inertia moment about A
I_A = ; % units
% Effective radius
r_eff=; 
% Effective mass
m_eff = ;
% Effective magnitude of the rotating force
F_0 = ;

% Time vector for theoretical response
t_theo = 0:1/fs:t_end-(1/fs);
phi = pi/2; % Starting phase

% apply predicted steady state solution 
x_steady_theo = ;

%% Plot experimental and theoretical data

% locate peaks
[peaks_exp,loc_exp] = findpeaks(x_for_res,'MinPeakDistance',(floor((1/4)*fs)-50),'MinPeakHeight',0.0005); %please do not edit this command

% trim experimental data so it starts at the same time as theoretical
t_exp_length = length(t_exp(loc_exp(1):end)); %please do not edit this command
t_exp_theo = linspace(0,(t_exp(end)-t_exp(loc_exp(1))),t_exp_length); %please do not edit this command

figure; % Students must complete this figure
grid on
% Insert plot variables for theoretical and experimental forced response here
% Insert x label here with proper units
% Insert y label here with proper units
% Insert legend here
