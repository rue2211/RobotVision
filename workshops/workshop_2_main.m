clc;
clear all;  % It's good practice to write clc and clear all at the start

%%
% READING CONSTANTS FILE TO LOAD CONSTANTS 
%     We will be using c, omega_ie, Omega_ie etc. from the
%     "Define_Constants.m" file
fid = fopen('Define_Constants.m');
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    eval(tline)
end

%% a) Initialise the Kalman filter state vector estimate & error cov
[x_est,P_matrix] = Initialise_GNSS_KF; %init state uncertainties 10m and clock offset is 0.1ms-1

%% b) Compute the transition matrix 
tmat = [eye(3) (1*eye(3)) zeros([3 1]) zeros([3 1]);
        zeros(3) eye(3) zeros([3 1]) zeros([3 1]);
        zeros([1 3]) zeros([1 3]) 1 1;
        zeros([1 3]) zeros([1 3]) 0 1];

%% c) Compute the system noise covariance matrix
sysNcov = [(1/3 * 5 * 1 * eye(3)) (1/2 * 5 * 1 * eye(3)) zeros([3 1]) zeros([3 1]);
            (1/2 * 5 * 1 * eye(3)) (5 * 1 * eye(3)) zeros([3 1]) zeros([3 1]);
            zeros([1 3]) zeros([1 3]) ((0.01 * 1) + (1/3 * 0.04 * 1)) (1/2 * 0.04 * 1);
            zeros([1 3]) zeros([1 3]) (1/2 * 0.04 * 1) (0.04 * 1)];

%% d) propogating state est
state_est = tmat * x_est;

%% e) propogate error cov mat
P_est = tmat * P_matrix * tmat' + sysNcov;

%% f) Predicting ranges and 
%  g) Computing los vector and
%  h) Predicting range rates

% rej --> position of satellite --> from satellite_position_and_velocity.m
% rea --> position estimate of receiver antenna --> x_est
% remember, this is for time stamp 0

sat_nos = [4 5	9	14	15	19	20	24	29	30]; % hard coding these for now
num_satellites = length(sat_nos);

sat_ECEF = zeros(6, num_satellites); % This is rej
for n = 1:num_satellites
    [pos, vel] = Satellite_position_and_velocity(0, sat_nos(n));
    sat_ECEF(:,n) = [transpose(pos); transpose(vel)];
end

receiver_ECEF = state_est; % This is rea

% Loop through each satellite
predicted_ranges = zeros(num_satellites, 1);
predicted_range_rates = zeros(num_satellites, 1);
u_los = zeros(3, num_satellites);

for n = 1:num_satellites
    % Initial range computation with CIe = identity matrix
    C_Ie = eye(3);
    temp = (C_Ie * sat_ECEF(1:3, n) - receiver_ECEF(1:3));
    initial_range = sqrt(temp'*temp); 
    %disp(['Predicted Range to Satellite ', num2str(sat_nos(n)), ' at Time 0 (Initial): ', num2str(initial_range)]);
    
    % Calculate the modified range with Sagnac correction
    C_Ie = [1, omega_ie * initial_range / c, 0; 
            -omega_ie * initial_range / c, 1, 0; 
            0, 0, 1];
    temp = (C_Ie * sat_ECEF(1:3, n) - receiver_ECEF(1:3));
    range_with_sagnac = sqrt(temp' * temp);
    %disp(['Predicted Range to Satellite ', num2str(sat_nos(n)), ' at Time 0 (Updated): ', num2str(range_with_sagnac)])
    
    % Store the predicted range
    predicted_ranges(n) = range_with_sagnac;
    
    % Line of sight vector
    u_los(:, n) = temp / predicted_ranges(n);
    
    % Predicting range rates
    temp2sat = sat_ECEF(4:6,n) + Omega_ie*sat_ECEF(1:3,n);
    temp2rec = receiver_ECEF(4:6) + Omega_ie*receiver_ECEF(1:3);
    predicted_range_rates = u_los'*(C_Ie*temp2sat - temp2rec);
    
end

%% i) Computing H matrix
H = zeros(num_satellites*2, 8);

for n = 1:num_satellites
    H(n,:) = [-u_los(1, n) -u_los(2, n) -u_los(3, n) 0 0 0 1 0];
end

for n=num_satellites+1:num_satellites*2
    i = n - num_satellites;
    H(n,:) = [0 0 0 -u_los(1, i) -u_los(2, i) -u_los(3, i) 0 1];
end

%% j) Compute measurement noise covariance matrix

% Given error standard deviations
std_dev_pseudorange = 10;  % meters
std_dev_pseudorange_rate = 0.05;  % m/s

% Measurement noise covariance matrix R
R = zeros(num_satellites * 2);

% diagonal elements for pseudorange measurements
R(1:num_satellites, 1:num_satellites) = eye(num_satellites) * std_dev_pseudorange^2;

% diagonal elements for pseudorange rate measurements
R(num_satellites+1:end, num_satellites+1:end) = eye(num_satellites) * std_dev_pseudorange_rate^2;

%% k) Kalman Gain
K = P_est * H' * inv(H * P_est * H' + R);

%% l) delta z
measured_pseudo_range = csvread('Workshop2_Pseudo_ranges.csv', 1, 1, [1 1 1 length(sat_nos)]);
measured_pseudo_range_rate = csvread('Workshop2_Pseudo_range_rates.csv', 1, 1, [1 1 1 length(sat_nos)]);

%x_minus = [receiver_ECEF; clock_offset];
delta_z_minus = zeros(num_satellites*2,1);

% Storing delta_zs for pseudo ranges (clock offset)
for n = 1:num_satellites
    delta_z_minus(n) = measured_pseudo_range(n) - predicted_ranges(n) - state_est(7);
end
% Storing delta_zs for pseudo range rates (clock drift)
for n = num_satellites+1:num_satellites*2
    i = n - num_satellites;
    delta_z_minus(n) = measured_pseudo_range_rate(i) - predicted_range_rates(i) - state_est(8);
end

%% m) updating state estimate
state_est = state_est + K*delta_z_minus;

%% n) updating error covariance
P_est = (eye(8) - K*H)* P_est;  %% This is perfectly correct

%% o) Conversion to ECEF
[latitude,longitude,height,vel]= pv_ECEF_to_NED(state_est(1:3), state_est(4:6));

% Converting to degrees
latitude = latitude* rad_to_deg;
longitude = longitude*rad_to_deg;

disp(['Latitude: ', num2str(latitude)]);
disp(['Longitude: ', num2str(longitude)]);
disp(['Height: ', num2str(height)]);
disp(['Velocity: ']);
vel

