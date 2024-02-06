function first_epoch_positioning_without_initialization()

    % Constants
    c = 299792458;  % Speed of light in m/s
    omega_ie = 7.2921159e-5;  % Earth's rotation rate in rad/s

    % Satellite numbers
    sat_nos = [2, 17, 18, 22, 23, 26, 27, 28];

    % Receiver initial position in cartesian ECEF coordinates
    receiver_ECEF = [0; 0; 0];

    % Receiver clock offset
    clock_offset = 1;

    % Load pseudo-ranges from CSV file
    measured_pseudo_range = csvread('Workshop1_Pseudo_ranges.csv', 1, 1, [1 1 1 length(sat_nos)]);

    % Call function to get satellite ECEF positions
    sat_ECEF = get_satellite_positions(sat_nos);

    % Iterative process for first-epoch positioning without initialization
    counter = 0;

    while true
        counter = counter + 1;

        % Call function to perform least squares computation
        [receiver_ECEF, clock_offset, converged] = least_squares_computation(receiver_ECEF, clock_offset, sat_ECEF, measured_pseudo_range, omega_ie, c);

        disp('Position and receiver clock offset are:')
        disp(num2str(receiver_ECEF'))
        disp(['Clock Offset: ', num2str(clock_offset)])

        if converged
            disp(['Total no. of iterations: ', num2str(counter)]);
            break;
        end

        disp(' ');

    end

end