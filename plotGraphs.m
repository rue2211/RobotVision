
% Define file names
file_names = {'GNSS_Solution_Output.csv', 'DR_Solution_Output.csv','Motion_Profile.csv'}; % Update with your file names

% Initialize plot
figure;
hold on;

% Loop through each file
for i = 1:length(file_names)
    % Load CSV file
    data = csvread(file_names{i}, 1); % Skip the first row
    
    % Extract data columns
    time = data(:, 1); % Assuming the first column contains time
    longitude = data(:, 3); % Assuming the second column contains longitude
    latitude = data(:, 2); % Assuming the third column contains latitude

    % Plot the data
    plot3(time, latitude, longitude);
end

% Customize plot
xlabel('Time');
ylabel('Longitude');
zlabel('Latitude');
title('Longitude and Latitude over Time');
grid on;
% Add legend if needed
legend('GNSS Only', 'DR Only', 'Integrated System'); % Update with appropriate legend entries

% Reset hold state
hold off;

