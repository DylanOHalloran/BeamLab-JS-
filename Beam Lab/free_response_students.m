        %%%%%%%%%%%%%%%% Theo&Exp Free response %%%%%%%%%%%%%%
        %%%%%%%%%%%% 3B5 Vibrating Beam Laboratory %%%%%%%%%%%
        %%%%%%% If in trouble, contact the demonstrator %%%%%%
        %%%%%%%%%% Federico Caruso (fcaruso@tcd.ie) %%%%%%%%%%

clc
clear
close all

fs = 1000;         % sample rate
%% Input data from the lab sheet
M = 1.633; % Concentrated mass, kg
m_cl_ds =0.4334 ; % Mass of clamps for spring and dashpot, kg
m_motor = 0.283; % Mass of the motor (including clamps and rotor), kg
m_beam = 1.94; % Mass of the beam, kg
l_1 = 0.34; % Length l1, m
l_2 = 0.482; % Length l2, m
l_3 = 0.505; % Length l3, m
l_4 = 0.62; % Length l4, m
l_beam = 0.76; % Length of the beam, m
k = 985; % Spring constant, N/m
I_beam = (1/3 * m_beam * l_beam ^ 2); % Inertia moment of the beam kg*m^2 (1/3 * m_beam * l_beam ^ 2)

% measurement factor
measurement_pt_factor = l_1/l_3; % please do not edit this command

%% Calculate and plot Experimental Response

% Load experimental data
exp_free_res = load('Free_Response_20th');

%% Separate into two vectors, change the measurement point to L1 and convert from mm to m
% Displacement
x_trans_exp = (exp_free_res(:,2) - mean(exp_free_res(1:500,2))) * 0.001;
x_trans_exp = x_trans_exp/l_3; % conversion from linear to rotational components
% Time
t_exp = exp_free_res(:,1);      %time variable separation


% plot the experimental data
figure; 
plot(t_exp,x_trans_exp,'r-')
xlabel('Time (s)'); ylabel('Displacement (mm)')
title('Raw Experimental Data')
grid on; 
hold on; 

% definition of the initial displacement
 [tempx,tempy] = ginput(1);
ind = find(t_exp>tempx);
ind = ind(1);
ave_disp = mean(x_trans_exp(ind:end));
free_exp_res = x_trans_exp - ave_disp;

plot(t_exp,free_exp_res,'b-')
grid on
xlabel('time (s)')
ylabel('Displacement (rad)')
hold off

%% locate peaks and ignore initial displacement (do not edit this block)
[peaks_exp,loc_exp] = findpeaks(free_exp_res,'MinPeakDistance',(floor((1/4)*fs)-50),'MinPeakHeight',0.001);
peaks_exp(1) = [];
loc_exp(1) = [];

%% fit a an exponential curve to peak data (do not edit this block)
delta_exp = fit(t_exp(loc_exp),peaks_exp,'exp1');
%Fourier Transform to determine the natural frequency of the system
[Pxx,F] = pwelch(free_exp_res(loc_exp(1):end),hanning(4000),[],length(free_exp_res),fs);

% plot data
figure;
subplot(311)
plot(t_exp,free_exp_res,t_exp(loc_exp),peaks_exp,'r+')
title('Fitting Experimental Data');
grid on
xlabel('Time (s)'); 
ylabel('Displacement (rad)')
subplot(312)
plot(delta_exp,'-b',t_exp(loc_exp),peaks_exp,'r+')
grid on
subplot(313)
semilogy(F,Pxx)
grid on
xlabel('Frequency (Hz)'); 
ylabel('Amplitude')
%% Natural frequency and damping ratio from experimental data
%% !! Students must determine damping ratio and natural frequency here !!
max_amp_th = max(Pxx);
f_n = 3.333 ; % Experimental natural frequency from Fast Fourier transform, Hz
%%^^ Above chosen as frequency at max amplitude ^^

w_n_exp = f_n * (2 * pi); % Natural frequency in rad/s
delta_coeff = coeffvalues(delta_exp); % Extract coefficient values

xi = -(delta_coeff(2) / w_n_exp) ; % Calculated Damping ratio from fitted data

%% Calculate Theoretical Response
% Theoretical Natural Frequency
k1 = (I_beam + M * (l_1 ^ 2) + m_motor * (l_3 ^ 2) + m_cl_ds * (l_4 ^ 2));
w_n_theo = sqrt((k * (l_4 ^ 2))/k1) ; % Theoretical natural Frequency

f_n_theo = w_n_theo*(1/(2*pi));

% inputs for theoretical model
B = peaks_exp(1);                                  % peak amplitude measured from experimental data
t_end = 6 ;                                         % test duration
t_theo = 0:1/fs:t_end-(1/fs);                      % time theoretical response

% calculate theoretical transient response
free_theo_res = zeros (1, t_end * 1000);



for i=1:length(free_theo_res)
free_theo_res(i)  = (B * exp(- xi * f_n_theo * t_theo(i)) *cos(t_theo(i) * f_n_theo * sqrt(1 - (xi ^ 2)))); % Theoretical free response of the beam
end
%locate peaks and fit an exponential curve to data
[peaks_theo,loc_theo] = findpeaks(free_theo_res,'MinPeakDistance',(floor((1/f_n)*fs)-50));
fit_theo = fit(t_theo(loc_theo)',peaks_theo','exp1');

% plot theoretical data
figure;
subplot(211)
plot(t_theo,free_theo_res,t_theo(loc_theo),peaks_theo,'r+')
title('Fitting Theoretical Data')
grid on
xlabel('Time (s)'); ylabel('Displacement (rad)')

subplot(212)
plot(fit_theo,'-b',t_theo(loc_theo),peaks_theo,'r+')
grid on

% determine theoretical damping ratio
xi_fit_theo = -fit_theo.b/(w_n_theo);

% Plot experimental and theoretical data

% trim experimental data so it starts at the same time as theoretical
t_exp_length = length(t_exp(loc_exp(1):end)); 
t_exp_theo = linspace(0,(t_exp(end)-t_exp(loc_exp(1))),t_exp_length); % time experimental response
figure; % Students must complete this figure
grid on
% Insert plot variables for theoretical and experimental free response here
plot(t_theo, free_theo_res, '-g', t_exp, free_exp_res, 'r-')
title('Theoretical and Experimental Free Response')

% Insert x label with appropriate units here
% Insert y label with appropriate units here
% Insert legend here