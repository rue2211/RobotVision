function results = W3Task1Results,('Workshop3_Speed_Heading.csv');

    h = 37.4; 
    
    %Latitude and longitude at time=0 in radians
    initial_latitude = deg2rad(50.4249580);
    initial_longitude = deg2rad(-3.5957974);

    %Load data file
    data = csvread('Workshop3_Speed_Heading.csv');

    %Initialise arrays for variables
    times = data(:, 1);
    latitude = zeros(size(data, 1), 1);
    longitude = zeros(size(data, 1), 1);
    v_N = zeros(size(data, 1), 1);
    v_E = zeros(size(data, 1), 1);

    latitude(1) = initial_latitude;
    longitude(1) = initial_longitude;

    for k = 2:size(data, 1)
        %Calculate average velocity
        avg_v_N = 0.5 * data(k, 2) * (cos(deg2rad(data(k, 3))) + cos(deg2rad(data(k-1, 3))));
        avg_v_E = 0.5 * data(k, 2) * (sin(deg2rad(data(k, 3))) + sin(deg2rad(data(k-1, 3))));

        %Calculate damped instantaneous DR velocity
        v_N(k) = 1.7 * avg_v_N - 0.7 * v_N(k-1);
        v_E(k) = 1.7 * avg_v_E - 0.7 * v_E(k-1);

        [R_N, R_E] = Radii_of_curvature(latitude(k-1));

        %Difference in time
        diff_time = times(k) - times(k-1);

        %Calculate latitude and longitude
        latitude(k) = latitude(k-1) + (v_N(k) * diff_time / (R_N + h));
        longitude(k) = longitude(k-1) + (v_E(k) * diff_time / ((R_E + h) * cos(latitude(k-1)))); 
    end

    %Convert latitude and longitude back to degrees
    latitude_deg = rad2deg(latitude);
    longitude_deg = rad2deg(longitude);

    %Create and display table for results
    results = table(times, latitude_deg, longitude_deg, v_N, v_E,...
        'VariableNames', {'Time', 'Latitude', 'Longitude', 'DR_Velocity_North', 'DR_Velocity_East'});

    disp(results);
end