
[x_est, P_matrix] = Initialise_GNSS_KF; % initialise with a set of zeros 
range = csvread("Pseudo_ranges.csv");
rangeepoch = range(2:end, 1);


% Initialize storage for results
num_epochs = length(rangeepoch);
% resultsTable = table('Size', [num_epochs, 6], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'Latitude', 'Longitude', 'Height', 'VelX', 'VelY', 'VelZ'});

% Initial state and covariance matrix
current_state = x_est;
current_P_est = P_matrix;

% Loop over epochs
for r = 1:num_epochs
    if r == num_epochs
        break;  % Exit the loop after the specified number of epochs
    end
    % Get the current pseudo-range measurement for the epoch
    current_pseudo_range = rangeepoch(r);

    % Call the Kalman filter function for the current epoch
    [current_state, current_P_est, latitude, longitude, height, vel] = GNSSImp(current_state, current_P_est, r);


end