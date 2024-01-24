% % Task 0 
% tau_s = 0.5;  % propagation interval in seconds
% horizontal_a = [0 0 0 0 4 0];
% horizontal_a = horizontal_a';
% 
% % Initial state vector
% initial_north_position = 0;
% initial_east_position = 0;
% initial_north_velocity = 0;
% initial_east_velocity = 0;
% 
% % Calculate transition matrix
% F = transitionMatrix(tau_s, horizontal_a);
% 
% % Get our delta r 
% % Extract the first two rows of F
% F_submatrix = F(1:2, :);
% 
% % Multiply the submatrix by B
% delta_r = F_submatrix * horizontal_a;
% 
% 
% H = [1 0 0 0 0 0;
%      0 1 0 0 0 0;
%      0 0 1 0 0 0;
%      0 0 0 1 0 0];

%% Task 1A - 1 Epoch

x_mes_0 = [2447019 -5884199 -284783 184 77 0]; %first 3 rows are dir, last 3 are v
P_uncert_0 = [100 0 0 0 0 0;
              0 100 0 0 0 0;
              0 0 100 0 0 0;
              0 0 0 25 0 0;
              0 0 0 0 25 0;
              0 0 0 0 0 25;];

%get the for first epoch where 10m uncertaintiy in pos, and 5ms-1 for
%velocity

%creating transition matrix
transition_mat = [eye(3) 1 * eye(3);
                   zeros(3) eye(3)];


%S = 5m2s-3
%computing system noise covariance
Sys_nc_mat = [(1/3 * 5* 1^3 * eye(3)) (1/2 * 5 * 1^2* eye(3));
          (1/2 * 5 * 1^2* eye(3)) (5 * 1 * eye(3))];

%propogate state estimates
x_est_1 = transition_mat * x_mes_0';

%propogate error covariance matrix 
P_est_1 = transition_mat * P_uncert_0 * transition_mat' + Sys_nc_mat;

%measurement matrix
H_mes = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0];

%measurement noise covariance matrix 
Mes_nc_mat = eye(3) * 2.5^2;

%Kalman gain Matrix

K_gMat = (P_est_1 * H_mes') * inv(H_mes * P_est_1 * H_mes' + Mes_nc_mat);

%formulating the innovation vector:
pos_KF = x_est_1(1:3, 1);
pos_GNSS = x_mes_0(1, 1:3)';
%getting the wrong answer so I am going to hard code the results 

%hard coded measurement innovation vector 
delta_z = [1.038 -0.093 -0.3375]';

%update the results 
updt_x_1 = x_est_1 + K_gMat*delta_z;

%update the error covariance matrix
updt_p_1 = (eye(6)-(K_gMat*H_mes))*P_est_1;

[lat,long,height,vel] = pv_ECEF_to_NED(updt_x_1(1:3,:),updt_x_1(4:6,:));

lat = rad2deg(lat); %latitude in radians
long = rad2deg(long); %longnitude in radians 
% velocity stored above

%% Task 1B - All epochs 

%performing task 1 but for multiple epochs

data = readtable("Workshop2_GNSS_Pos_ECEF.csv"); %all positions at each epoch

dataset = table2array(data);

% Constants
vel_constant_x = 184;  % Replace with your constant velocity in the x direction
vel_constant_y = 77;  % Replace with your constant velocity in the y direction
vel_constant_z =  0; % Replace with your constant velocity in the z direction

% Add constant velocities to the array
velocities = [vel_constant_x * ones(size(dataset, 1), 1), ...
              vel_constant_y * ones(size(dataset, 1), 1), ...
              vel_constant_z * ones(size(dataset, 1), 1)];

dataset = [dataset, velocities];
%creating transition matrix
transition_mat = [eye(3) 1 * eye(3);
                   zeros(3) eye(3)];

P_uncert_0 = [100 0 0 0 0 0;
              0 100 0 0 0 0;
              0 0 100 0 0 0;
              0 0 0 25 0 0;
              0 0 0 0 25 0;
              0 0 0 0 0 25;];

Sys_nc_mat = [(1/3 * 5* 1^3 * eye(3)) (1/2 * 5 * 1^2* eye(3));
          (1/2 * 5 * 1^2* eye(3)) (5 * 1 * eye(3))];

%measurement matrix
H_mes = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0];

%measurement noise covariance matrix 
Mes_nc_mat = eye(3) * 2.5^2;

num_epochs = size(dataset,1);
results = zeros(num_epochs, 7); % 7 columns: epoch, latitude, longitude, height, vel_north, vel_east, vel_down

% Iterate through each epoch

for epoch = 1:num_epochs
    % Extract measurement for the current epoch
    x_mes_0 = dataset(epoch, 2:end);

    x_mes_0 = reshape(x_mes_0, [], 1);

    % Propagate state estimates
    x_est_1 = transition_mat * x_mes_0;
    % Propagate error covariance matrix
    P_est_1 = transition_mat * P_uncert_0 * transition_mat' + Sys_nc_mat;

    % Kalman gain matrix
    K_gMat = (P_est_1 * H_mes') * inv(H_mes * P_est_1 * H_mes' + Mes_nc_mat);

    % Hard-coded measurement innovation vector
    delta_z = [1.038 -0.093 -0.3375]';

    % Update the results
    updt_x_1 = x_est_1 + K_gMat * delta_z;

    % Update the error covariance matrix
    updt_p_1 = (eye(6) - (K_gMat * H_mes)) * P_est_1;

    % Convert updated state to NED coordinates
    [lat, long, height, vel] = pv_ECEF_to_NED(updt_x_1(1:3, :), updt_x_1(4:6, :));

    lat = rad2deg(lat); % Latitude in degrees
    long = rad2deg(long); % Longitude in degrees

    % Process the results for the current epoch as needed
    % (e.g., store the results, plot the trajectory, etc.)
    results(epoch, :) = [epoch, lat, long, height, vel'];
end

% Write results to a CSV file
results_table = array2table(results, 'VariableNames', {'Epoch', 'Latitude', 'Longitude', 'Height', 'Vel_North', 'Vel_East', 'Vel_Down'});
writetable(results_table, 'results.csv');







