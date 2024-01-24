

[x_est,P_matrix] = Initialise_GNSS_KF; %init state uncertainties 10m and clock offset is 0.1ms-1

%transition matrix

tmat = [eye(3) (1*eye(3)) zeros([3 1]) zeros([3 1]);
        zeros(3) eye(3) zeros([3 1]) zeros([3 1]);
        zeros([1 3]) zeros([1 3]) 1 1;
        zeros([1 3]) zeros([1 3]) 0 1];

%system noise covariance

sysNcov = [(1/3 * 5 * 1 * eye(3)) (1/2 * 5 * 1 * eye(3)) zeros([3 1]) zeros([3 1]);
            (1/2 * 5 * 1 * eye(3)) (5 * 1 * eye(3)) zeros([3 1]) zeros([3 1]);
            zeros([1 3]) zeros([1 3]) ((0.01 * 1) + (1/3 * 0.04 * 1)) (1/2 * 0.04 * 1);
            zeros([1 3]) zeros([1 3]) (1/2 * 0.04 * 1) (0.04 * 1)];

%propogating state est

state_est = tmat * x_est;

%propogate error cov mat

cov_mat = tmat * P_matrix * tmat' + sysNcov;

%range prediction - via sagnac 
%constants
c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s

%sagnac w identity - before filling out

SEC_mat_I = eye(3);

%range w identity, need to call the pseudo ranges like in exercises before
%need a way to concatinate the ranges 
M = readmatrix("Workshop2_Pseudo_ranges.csv");
v = readmatrix("Workshop2_Pseudo_range_rates.csv");
%extract satelite numbers
flightnum = M(1,2:10);


% Initialize a cell array to store results
results = cell(numel(flightnum), 4);


% Loop through the columns
for idx = 1:numel(flightnum)
    % Extract satelite number and data for the current column
    flightNumber = flightnum(idx);

    % Call the function with the current value of j_num and t
    outb = Satellite_position_and_velocity(v, flightNumber);
    
    % Store flight number and corresponding result in the cell array
    results{idx, 1} = flightNumber;
        % Use a loop to assign each element of outb individually
    for i = 1:numel(outb)
        results{idx, i + 1} = outb(i);
    end
end

%rej -> pseudo ranges
%rea -> maybe pos for state estimate
%need to loop over and find the differences between each satelite to the
%user position 
rea = state_est(1:3, 1);

% Initialize an array to store differences
% Loop through satellites and calculate differences
r_diff = cell(1, numel(flightnum));
raj_id = cell(1, numel(flightnum));
Sag_raj = cell(1, numel(flightnum));
raj_acc = cell(1, numel(flightnum));
u_aj = cell(1, numel(flightnum));

for sat_idx = 1:numel(flightnum)
    % Extract satellite position for all epochs
    rej = cell2mat(results(sat_idx, 2:end)');

    % Calculate the differences between satellite position and user position
    differences = rej - rea;

    % Store the result for the current satellite
    r_diff{sat_idx} = differences;
    
    raj_id{sat_idx} = sqrt((SEC_mat_I * differences)' * (SEC_mat_I * differences));

    % Calculate Sag_raj matrix for the current satellite
    Sag_raj{sat_idx} = [1, (omega_ie * raj_id{sat_idx}/c), 0;
                        - (omega_ie * raj_id{sat_idx}/c), 1, 0;
                        0, 0, 1];

    raj_acc{sat_idx} = sqrt((Sag_raj{sat_idx} * differences)' * (Sag_raj{sat_idx} * differences));

    u_aj{sat_idx} = (Sag_raj{sat_idx} * differences)/raj_acc{sat_idx};
end





