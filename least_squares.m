function [receiver_ECEF, clock_offset, converged] = least_squares(receiver_ECEF, clock_offset, sat_ECEF, measured_pseudo_range)
    % Function to perform least squares computation
    
    c = 299792458; % Speed of light in m/s
    omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
    converged = false;

    while(true)
    % Initialize variables
    predicted_ranges = zeros(1, length(sat_ECEF));
    u_los = zeros(3, length(sat_ECEF));
    delta_z_minus = zeros(length(sat_ECEF), 1);
    H = zeros(length(sat_ECEF), 4);

    for n = 1:length(sat_ECEF)

        % Initial range computation with CIe = identity matrix
        C_Ie = eye(3);
        temp = (C_Ie * sat_ECEF(:, n)) - receiver_ECEF;
        initial_range = sqrt(temp' * temp);

        % Final range computation with the updated CIe using the range predicted above
        C_Ie = [1 omega_ie * initial_range / c 0; -omega_ie * initial_range / c 1 0; 0 0 1];
        temp = (C_Ie * sat_ECEF(:, n)) - receiver_ECEF;
        range_with_sagnac = sqrt(temp' * temp);

        predicted_ranges(n) = range_with_sagnac;
        u_los(:, n) = temp / predicted_ranges(n);

        delta_z_minus(n) = measured_pseudo_range(n) - predicted_ranges(n) - clock_offset;
        H(n, :) = [-u_los(1, n) -u_los(2, n) -u_los(3, n) 1];

    end

    % Compute updated receiver ECEF and clock offset
    updated_receiver_ECEF = [receiver_ECEF; clock_offset] + inv(H' * H) * H' * delta_z_minus;

    if norm(receiver_ECEF - updated_receiver_ECEF(1:3)) <= 0.1
        converged = true;
        break;
    end

    receiver_ECEF = updated_receiver_ECEF(1:3);
    clock_offset = updated_receiver_ECEF(4);

end