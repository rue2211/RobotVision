function sat_ECEF = get_satellite_positions(sat_nos)
    % Function to compute ECEF positions of satellites

    sat_ECEF = zeros(3, length(sat_nos));

    for n = 1:length(sat_nos)
        sat_ECEF(:, n) = transpose(Satellite_position_and_velocity(0, sat_nos(n)));
    end

end
