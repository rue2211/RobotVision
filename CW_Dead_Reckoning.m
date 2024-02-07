clc;
clear all;
Define_Constants;

DR_inputs = readmatrix('Dead_reckoning.csv');

% time interval
t =0.5; 

% Initial h, latitude and longitude at time 0
h = 39.2043;
GNSS_lat_b = 51.5092543897043*deg_to_rad;
GNSS_long_b = -0.161045151548226*deg_to_rad;

% Correcting Speed Sensor values using known errors
speed_scale_factor = 0.03; 
speed_noise_std = 0.05;
speed_quant = 0.02;
DR_vel_c = zeros(851,1); % Corrected velocities using the rear(driver) wheels
for k = 1:851
    if DR_inputs(k,4) >= 0 && DR_inputs(k,5) >= 0
        DR_vel_c(k,1) = (DR_inputs(k,4) - (speed_noise_std + speed_quant + DR_inputs(k,4)*speed_scale_factor) + DR_inputs(k,5) - (speed_noise_std + speed_quant + DR_inputs(k,5)*speed_scale_factor))/2;

    elseif DR_inputs(k,4) >= 0 && DR_inputs(k,5) < 0
            DR_vel_c(k,1) = (DR_inputs(k,4) - (speed_noise_std + speed_quant + DR_inputs(k,4)*speed_scale_factor) + DR_inputs(k,5) + (speed_noise_std + speed_quant - DR_inputs(k,5)*speed_scale_factor))/2;
    
    elseif DR_inputs(k,4) < 0 && DR_inputs(k,5) < 0
            DR_vel_c(k,1) = (DR_inputs(k,4) + (speed_noise_std + speed_quant - DR_inputs(k,4)*speed_scale_factor) + DR_inputs(k,5) + (speed_noise_std + speed_quant - DR_inputs(k,5)*speed_scale_factor))/2;
    else
        DR_vel_c(k,1) = (DR_inputs(k,4) + (speed_noise_std + speed_quant - DR_inputs(k,4)*speed_scale_factor) + DR_inputs(k,5) - (speed_noise_std + speed_quant + DR_inputs(k,5)*speed_scale_factor))/2;
    end

 end


 % Initial velocities at time 0
 DR_vel_N = DR_vel_c(1,1)*cos(DR_inputs(1,7)*deg_to_rad);
 DR_vel_E = DR_vel_c(1,1)*sin(DR_inputs(1,7)*deg_to_rad);
 DR_vel_NE_c = [DR_vel_N + (speed_noise_std + speed_quant - DR_vel_N * speed_scale_factor);
             DR_vel_E - (speed_noise_std + speed_quant + DR_vel_E * speed_scale_factor)];

 GNSS_lat = GNSS_lat_b;
 GNSS_long = GNSS_long_b;
 heading = DR_inputs(1,7)*deg_to_rad;

 DR_result = zeros(851,5);
 DR_result(1,1) = 0;
 DR_result(1,2) = GNSS_lat_b*rad_to_deg;
 DR_result(1,3) = GNSS_long_b*rad_to_deg;
 DR_result(1,4) = DR_vel_NE_c(1);
 DR_result(1,5) = DR_vel_NE_c(2);

 % Constants for gyroscope measurement correction
 gyr_bias_std = 1;
 gyr_scale_factor = 0.01 + 0.001; %scale factor and cross coupling error std
 gyr_quant = 0.0002;


 for k = 1:850

 % Correcting the gyroscope values using known errors
 if DR_inputs(k+1,6)>=0
     omega = DR_inputs(k+1,6) - (gyr_bias_std*deg_to_rad + gyr_quant + DR_inputs(k+1,6)*gyr_scale_factor);
 else
     omega = DR_inputs(k+1,6) + (gyr_bias_std*deg_to_rad + gyr_quant - DR_inputs(k+1,6)*gyr_scale_factor);
 end

 % Heading weight (smoothened heading)
 w = (1 + (0.01+0.001) * DR_inputs(k+1,6)*rad_to_deg + 0.0002*rad_to_deg)*0.5/4;
 heading_c = w*DR_inputs(k+1,7)*deg_to_rad +(1-w)*(heading+t*omega);

 %compute the overall speed
 avg_v_N = 0.5*(cos(heading_c)+cos(heading))*DR_vel_c(k,1);
 avg_v_E = 0.5*(sin(heading_c)+sin(heading))*DR_vel_c(k,1);
 v_ned = [avg_v_N; avg_v_E];

 [R_N,R_E]= Radii_of_curvature(GNSS_lat);
 %longitude and latitude
 GNSS_lat_c = GNSS_lat+(avg_v_N*t)/(R_N+h);
 GNSS_long_c = GNSS_long+(avg_v_E*t)/((R_E+h)*cos(GNSS_lat_c));


 %instantaneous velocity
 vel_n_and_e =1.7*[avg_v_N; avg_v_E]-0.7*DR_vel_NE_c; % Note that this is column vector [a;b] of size 2x1

 DR_result(k+1,1)= k*t ; %time
 DR_result(k+1,2)=GNSS_lat_c*rad_to_deg; %latitude
 DR_result(k+1,3)=GNSS_long_c*rad_to_deg; %longitude
 DR_result(k+1,4:5)=vel_n_and_e.'; %velocity in North and East
 DR_result(k+1,6)= heading_c; %Heading
 GNSS_lat =GNSS_lat_c;
 GNSS_long =GNSS_long_c;
 DR_vel_NE_c =vel_n_and_e;
 heading =heading_c;
 end

writematrix(DR_result,'DR_solution_Output.csv');
