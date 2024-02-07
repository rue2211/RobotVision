clc;
clear all;
Define_Constants;

%% 

% Obtaining DR Solution 
DR_soln = readmatrix("DR_solution_Output.csv");
% Obtaining GNSS Measurements 
GNSS_soln = readmatrix('GNSS_Solution_Output.csv');
% Initialing result matrix
final_output = zeros(851,5); % replace size here

% Storing columns into variables for readability in the operations within
% the KF loop
GNSS_lat = GNSS_soln(:,2)*deg_to_rad;
GNSS_long = GNSS_soln(:,3)*deg_to_rad;
GNSS_vel_N = GNSS_soln(:,5);
GNSS_vel_E = GNSS_soln(:,6);
h = GNSS_soln(:,4);

DR_lat = DR_soln(:,2)*deg_to_rad;
DR_long = DR_soln(:,3)*deg_to_rad;
DR_vel_N = DR_soln(:,4);
DR_vel_E = DR_soln(:,5);
DR_heading = DR_soln(:,6);

%%

% Initial error covariance matrix P
GNSS_lat_b = GNSS_lat(1);
GNSS_long_b = GNSS_long(1);
h_b = h(1);
[R_N,R_E]= Radii_of_curvature(GNSS_lat_b);

sig_v = 0.1; % error/uncertainty in velocity 
sig_r = 10; % error/uncertainty in position 
P = [sig_v^2   0     0  0;
     0      sig_v^2  0  0;
     0      0      sig_r^2/(R_N + h_b)^2  0;
     0     0       0    sig_r^2/((R_E + h_b)^2*(cos(GNSS_lat_b))^2)];

% Initial state_vec
state_vec = zeros(4,1);

% Known errors 
sdr = 0.2; % for noise covariance matrix (noise due to acceleration)
t = 0.5; % time difference tau
sig_gr = 5;% position measuremnet error standard deviation
sig_gv = 0.02;% velocity measuremnet error standard deviation


%% 
% KF loop

for i = 1:849
    
    %%% ------- PROPAGATION STEP ------%%%
    GNSS_lat_b = GNSS_lat(i);
    GNSS_long_b = GNSS_long(i);
    h_b = h(i);
    [R_N,R_E]= Radii_of_curvature(GNSS_lat_b);

    %transition matrix
    T = [1 0 0 0;
         0 1 0 0;
         0.5/(R_N + h_b) 0 1 0;
         0 0.5/((R_E + h_b)*cos(GNSS_lat_b)) 0 1];

    %noise covariance matrix
    Q = [sdr*t  0  0.5*(sdr*t^2)/(R_N + h_b) 0;
          0    sdr*t    0   0.5 * ((sdr)*t^2) / ( (R_E + h_b) * cos(GNSS_lat_b) )
          0.5*(sdr*t^2)/(R_N + h_b)  0  1/3*(sdr*t^3)/(R_N + h_b)^2  0;
          0   0.5 * ((sdr)*t^2) / ( (R_E + h_b) * cos(GNSS_lat_b))   0  1/3 * ((sdr)*t^3) / ( (R_E + h_b)^2 * cos(GNSS_lat_b)^2) ];

    % propagation step
    propagated_state_vec = T*state_vec;
    propagated_P = T*P*T.' + Q;

    %%% --------------------------------%%%



    %%% ------- MEASUREMENT STEP ------%%%
    % MEasurement matrix
    H = [0 0 -1 0;
        0 0 0 -1;
        -1 0 0 0;
        0 -1 0 0];
    
    % R - Measurement Noise Covariance Matrix
    GNSS_lat_c = GNSS_lat(i);
    GNSS_long_c = GNSS_long(i);
    h_c = h(i);

    [R_N,R_E]= Radii_of_curvature(GNSS_lat_c);
    R = [((sig_gr^2)/(R_N+h_c)^2) 0 0 0;
         0 ((sig_gr^2)/(R_E+h_c)^2*cos(GNSS_lat_c)^2) 0 0;
         0 0 sig_gv^2 0;
         0 0 0 sig_gv^2;];

    % Kalman gain matrix
    K = propagated_P*H.'*inv(H*propagated_P*H.'+R);

    % Measurement innovation vector
   
    del_z_inter = [GNSS_lat_c - DR_lat(i);
                   GNSS_long_c - DR_long(i);
                   GNSS_vel_N(i) - DR_vel_N(i);
                   GNSS_vel_E(i) - DR_vel_E(i)];

    del_z = del_z_inter - H*propagated_state_vec;
   

    % Update the state estimation and P matrix
    state_vec = propagated_state_vec + K*del_z;
    P = (eye(4) - K*H)*propagated_P;

    x_to_store = [state_vec(1);
                  state_vec(2);
                  state_vec(3)*rad_to_deg;
                  state_vec(4)*rad_to_deg];

    final_output(i,1) = (i-1)*t;
    final_output(i,2:5) = DR_soln(i,2:5).' + H*x_to_store;
     %%% --------------------------------%%%
end

integrated_GNSS_DR_output = [final_output DR_heading]
writematrix(integrated_GNSS_DR_output,'Motion_Profile.csv');

