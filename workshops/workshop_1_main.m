clc;
clear all;

%%
% READING CONSTANTS FILE TO LOAD CONSTANTS
fid = fopen('Define_Constants.m');
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    eval(tline)
end

%% 1A) First-epoch Positioning with Initialisation 

%% 1A) 1. To compute the ecef position of receiver antenna

% INPUTS
latitude = -33.821075;
longitude = 151.188496;
height = 120;
velocity = 0.0;

lat_r = latitude * deg_to_rad;
long_r = longitude * deg_to_rad;

receiver_ECEF = pv_NED_to_ECEF(lat_r, long_r, height, velocity);

%% 1A) 2. To compute the ecef position of satellite

sat_nos = [2, 17, 18, 22, 23, 26, 27, 28];
sat_ECEF = zeros(3, length(sat_nos));

time = 0;

for n = 1:length(sat_nos)
    sat_ECEF(:, n) = transpose(Satellite_position_and_velocity(0, sat_nos(n)));
end
%% 1A) 3. To predict the range
%      4. Compute line of sight vector

predicted_ranges = zeros(1, length(sat_nos));
u_los = zeros(3, length(sat_nos));

for n = 1:length(sat_nos)
    
    % Initial range computation with CIe = identity matrix
    C_Ie = eye(3);
    temp = (C_Ie * sat_ECEF(:, n)) - receiver_ECEF;
    initial_range = sqrt(temp'*temp); 
    disp(['Predicted Range to Satellite ', num2str(sat_nos(n)), ' at Time 0 (Initial): ', num2str(initial_range)]);
    
    % Final range computation with the updated CIe using the range
    % predicted above
    C_Ie = [1 omega_ie*initial_range/c 0; -omega_ie*initial_range/c 1 0; 0 0 1];
    temp = (C_Ie * sat_ECEF(:, n)) - receiver_ECEF;
    range_with_sagnac = sqrt(temp'*temp); 
    disp(['Predicted Range to Satellite ', num2str(sat_nos(n)), ' at Time 0 (Updated): ', num2str(range_with_sagnac)]);

    predicted_ranges(n) = range_with_sagnac;
    u_los(:, n) = temp / predicted_ranges(n);
    


end


%% 1A) 5. Formulating x-, H and delta z-

clock_offset = 1;

measured_pseudo_range = csvread('Workshop1_Pseudo_ranges.csv', 1, 1, [1 1 1 length(sat_nos)]);

x_minus = [receiver_ECEF; clock_offset];
delta_z_minus = zeros(length(sat_nos),1);
H = zeros(length(sat_nos), 4);


for n = 1:length(sat_nos)
    delta_z_minus(n) = measured_pseudo_range(n) - predicted_ranges(n) - clock_offset;
    H(n,:) = [-u_los(1, n) -u_los(2, n) -u_los(3, n) 1];
end


%% 1A) 6. Computing for x+ using unweighted least squares  %nsolve

x_plus = x_minus + inv(H'*H)*H'*delta_z_minus;
disp('Position and receiver clock offset are:')
for i=1:4
   disp(num2str(x_plus(i)))
end


%% 1A) 7. Computing for x+ using unweighted least squares  %nsolve

[latitude,longitude,height,vel]= pv_ECEF_to_NED(x_plus(1:3),0);

% Converting to degrees
latitude = latitude* rad_to_deg;
longitude = longitude*rad_to_deg;

disp(['Latitude: ', num2str(latitude)]);
disp(['Longitude: ', num2str(longitude)]);
disp(['Height: ', num2str(height)]);

%% 1B) First-epoch Positioning without Initialisation 

receiver_ECEF = [0;0;0];
counter = 0;

while(true)
    counter = counter + 1;
    x_minus = [receiver_ECEF; clock_offset];
    delta_z_minus = zeros(length(sat_nos),1);
    H = zeros(length(sat_nos), 4);
    
    predicted_ranges = zeros(1, length(sat_nos));
    u_los = zeros(3, length(sat_nos));
    
    for n = 1:length(sat_nos)
        
        % Initial range computation with CIe = identity matrix
        C_Ie = eye(3);
        temp = (C_Ie * sat_ECEF(:, n)) - receiver_ECEF;
        initial_range = sqrt(temp'*temp);     
        % Final range computation with the updated CIe using the range
        % predicted above
        C_Ie = [1 omega_ie*initial_range/c 0; -omega_ie*initial_range/c 1 0; 0 0 1];
        temp = (C_Ie * sat_ECEF(:, n)) - receiver_ECEF;
        range_with_sagnac = sqrt(temp'*temp); 
    
        predicted_ranges(n) = range_with_sagnac;
        u_los(:, n) = temp / predicted_ranges(n);
    
        delta_z_minus(n) = measured_pseudo_range(n) - predicted_ranges(n) - clock_offset;
        H(n,:) = [-u_los(1, n) -u_los(2, n) -u_los(3, n) 1];
    
    end
    
    updated_receiver_ECEF = x_minus + inv(H'*H)*H'*delta_z_minus;
    
    disp('Position and receiver clock offset are:')
    for i=1:4
       disp(num2str(updated_receiver_ECEF(i)))
    end
    
    
    if norm(receiver_ECEF - updated_receiver_ECEF(1:3))<= 0.1
        disp(" ");
        disp("Total no. of iterations " + counter);
        break;
    else
        receiver_ECEF = updated_receiver_ECEF(1:3);
    end

disp(' ');

end


%% 1C) First-epoch Positioning without Initialisation 
