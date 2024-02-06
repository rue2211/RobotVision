range = csvread('Pseudo_ranges.csv');



sat_nos = range(1, 2:end);

x_0 = zeros(8,1);
P_0 = diag([1^2, 1^2, 1^2, 0.1^2, 0.1^2 0.1^2, 100000^2, 200^2]);

sat_ecef = get_satellite_positions(sat_nos);
    
measured_pseudo_range_row = range(2, 2:end);
    
[x_f, clock_off] = least_squares(x_0(1:3), x_0(7), sat_ecef, measured_pseudo_range_row);


x_0(1:3) = x_f;
x_0(7) = clock_off;

range_epoch = range(2:end, 1);

num_epochs = length(range_epoch);

current_state = x_0;
current_P_est = P_0;

%[state_est1, P_est1, latitude, longitude, height, vel] = GNSS_KF(current_state, current_P_est, 2);

% Initialize arrays to store results
latitude_array = zeros(num_epochs, 1);
longitude_array = zeros(num_epochs, 1);
height_array = zeros(num_epochs, 1);
velocity_x_array = zeros(num_epochs, 1);
velocity_y_array = zeros(num_epochs, 1);
velocity_z_array = zeros(num_epochs, 1);

for r = 1:num_epochs 
    if r == 1 || r ==2
        continue;  % Exit the loop after the specified number of epochs
    end
    [x_r, clock_off_r] = least_squares(current_state(1:3), current_state(7), sat_ecef, measured_pseudo_range_row);
    % Call the Kalman filter function for the current epoch
    current_state(1:3) = x_r;
    current_state(7) = clock_off_r;
    [state_est1, P_est1, latitude, longitude, height, velocity] = GNSS_KF(current_state, current_P_est, r);
    % Store results in arrays
    latitude_array(r) = latitude;
    longitude_array(r) = longitude;
    height_array(r) = height;
    velocity_x_array(r) = velocity(1);
    velocity_y_array(r) = velocity(2);
    velocity_z_array(r) = velocity(3);
    % Update current state and covariance matrix
    current_state = state_est1;
    current_P_est = P_est1;
end

% Create a results table
results_table = table(range_epoch(3:end), latitude_array(3:end), longitude_array(3:end), height_array(3:end), velocity_x_array(3:end), velocity_y_array(3:end), velocity_z_array(3:end), 'VariableNames', {'Time', 'Latitude', 'Longitude', 'Height', 'Velocity_North', 'Velocity_East', 'Velocity_Down'});
disp(results_table);


