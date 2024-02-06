function [state_est1, P_est1, latitude, longitude, height, vel] = GNSS_KF(prev_state, prev_P_est, r)
%% loading information

%pseudo ranges and rates
pseudo_range = csvread('Pseudo_ranges.csv');
pseudo_range_rate= csvread('Pseudo_range_rates.csv');

%taking the values without the satellite numbers
pr = pseudo_range(r, 2:end); % r needs to start at row 2
prr = pseudo_range_rate(r, 2:end); 


deg_to_rad = 0.01745329252; % Degrees to radians conversion factor
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor
c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]);

    %% a) Initialise the Kalman filter state vector estimate & error cov
    % Previous state and covariance matrix
    if nargin < 1
        % If not provided, initialize with zeros or default values
        prev_state = zeros(8, 1);
    end

    if nargin < 2
        % If not provided, initialize with identity matrix or default values
        prev_P_est = eye(8);
    end

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
    state_est = tmat * prev_state;

    %% e) propogate error cov mat
    P_est = tmat * prev_P_est * tmat' + sysNcov;

    sat_nos = pseudo_range(1, 2:end);
    num_satellites = length(sat_nos);
    
    % Initialize storage for satellite positions, velocities, and receiver position
    sat_ECEF = zeros(6, num_satellites);
    receiver_ECEF = state_est; % This is rea, assuming state_est is defined earlier
    
    % Initialize storage for predicted ranges, range rates, and LOS vectors
    predicted_ranges = zeros(num_satellites, 1);
    predicted_range_rates = zeros(num_satellites, 1);
    u_los = zeros(3, num_satellites);
    
    % Compute satellite positions, velocities, and receiver's ECEF position
        for n = 1:num_satellites
        [pos, vel] = Satellite_position_and_velocity(r, sat_nos(n));
        sat_ECEF(:, n) = [transpose(pos); transpose(vel)];
        end
        
        % Compute the line-of-sight unit vector and predicted range rates
        for n = 1:num_satellites
        % Calculate the line-of-sight vector using the approximation method
        r_ej = sat_ECEF(1:3, n);
        r_ea = receiver_ECEF(1:3);
        
        % Initial range computation without Earth's rotation effect
        initial_range = norm(r_ej - r_ea);
        
        % Sagnac effect correction matrix
        C_Ie = [1, omega_ie * initial_range / c, 0; 
        -omega_ie * initial_range / c, 1, 0; 
        0, 0, 1];
        
        % Correct the satellite position with the Sagnac effect
        r_ej_corrected = C_Ie * r_ej;
        
        % Compute the line-of-sight vector
        u_los(:, n) = (r_ej_corrected - r_ea) / norm(r_ej_corrected - r_ea);
        
        % Store the predicted range
        predicted_ranges(n) = norm(r_ej_corrected - r_ea);
        
        % Predicting range rates
        v_ej = sat_ECEF(4:6, n); % Satellite velocity
        v_ea = receiver_ECEF(4:6); % Receiver velocity
        
        % Apply Sagnac correction to satellite velocity
        v_ej_corrected = C_Ie * (v_ej + Omega_ie * r_ej);
        
        % Predict the range rate
        predicted_range_rates(n) = u_los(:, n)' * (v_ej_corrected - (v_ea + Omega_ie * r_ea));
        end
    
    % The predicted_ranges and predicted_range_rates now contain the values for each satellite
    
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
    
    mat_sigp = eye(num_satellites) * 10^2;
    mat_sigr = eye(num_satellites) * 0.05^2;
    
    R = [mat_sigp zeros(num_satellites,num_satellites);
    zeros(num_satellites,num_satellites) mat_sigr];
    
    %% k) Kalman Gain
    K = P_est * H' * inv(H * P_est * H' + R);
    
    %% l) delta z

    % measured_pseudo_range_row = csvread('Pseudo_ranges.csv', r, 1, [r, 1, r, length(sat_nos)]);
    % measured_pseudo_range_rate_row = csvread('Pseudo_range_rates.csv', r, 1, [r, 1, r, length(sat_nos)]);
    
    %x_minus = [receiver_ECEF; clock_offset];
    delta_z_minus = zeros(num_satellites*2,1);
    
    % Storing delta_zs for pseudo ranges (clock offset)
    for n = 1:num_satellites
    % Assuming sat_nos corresponds to columns in the order they appear in the CSV,
    % and measured_pseudo_range_row directly gives the range values for those satellites.
    delta_z_minus(n) = pr(n) - predicted_ranges(n) - state_est(7);
    end

    % Storing delta_zs for pseudo range rates (clock drift)
    for n = num_satellites+1:num_satellites*2
        i = n - num_satellites;
        delta_z_minus(n) = prr(i) - predicted_range_rates(i) - state_est(8);
    end

    % Output the results
    state_est1 = state_est + K * delta_z_minus;
    P_est1 = (eye(8) - K * H) * P_est;

    % Conversion to ECEF
    [latitude, longitude, height, vel] = pv_ECEF_to_NED(state_est(1:3), state_est(4:6));

    % Converting to degrees
    latitude = latitude * rad_to_deg;
    longitude = longitude * rad_to_deg;

    
    % disp(['Epoch: ', num2str(r)]);
    % disp(['Latitude: ', num2str(latitude)]);
    % disp(['Longitude: ', num2str(longitude)]);
    % disp(['Height: ', num2str(height)]);
    % disp(['Velocity: ']);
    % vel

end
