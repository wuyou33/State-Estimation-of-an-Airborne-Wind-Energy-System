function [measurements, TA, time_s] = readArduinoData(path, configFilename, dataFilename, g_mps2)
plotting = 0; % 1: plotting delta_t

disable_magnetometer = 0;
if disable_magnetometer
    warning('available magnetometer is set to zero for all values');
end

run(fullfile(path, configFilename));

M = csvread(fullfile(path, dataFilename));
M = M'; % columnvector

[rowcount, ~] = size(M');

poti = find(contains(columns, 'potentiometer'));
measurements.meas.poti = M(poti, :);
% use first value to determine offset!
if length(measurements.meas.poti) > 0
    measurements.meas.poti = measurements.meas.poti - measurements.meas.poti(1);
end
minimum = min(measurements.meas.poti); % rotated clockwise by 90°, seen from above
% assumption: resistor value is proportial to the rotated angle
measurements.meas.poti_deg = measurements.meas.poti * 90/minimum;

x = find(contains(columns, 'acX'));
y = find(contains(columns, 'acY'));
z = find(contains(columns, 'acZ'));
measurements.meas.a_b = M([x, y, z], :)/g_mps2 * -1; % [m/s^2] --> [g], don't know why the direction is wrong
measurements.available.acc = ones(1, rowcount);

x = find(contains(columns, 'gyX'));
y = find(contains(columns, 'gyY'));
z = find(contains(columns, 'gyZ'));
measurements.meas.w_b = M([x, y, z], :);
display('Change gyro unit from °/s to rad/s');
measurements.meas.w_b = measurements.meas.w_b *pi/180;
measurements.available.gyr = ones(1, rowcount);

x = find(contains(columns, 'magX'));
y = find(contains(columns, 'magY'));
z = find(contains(columns, 'magZ'));
mag_temp = M([x, y, z], :);
measurements.available.mag = ones(1, rowcount);
mag_norm = vecnorm(mag_temp);
measurements.available.mag(mag_norm == 0) = 0;
%measurements.meas.m_b = mag_temp;
measurements.meas.m_b = mag_temp./mag_norm; % normalize
measurements.meas.m_b_unnormalized = mag_temp;
if disable_magnetometer
    measurements.available.mag = zeros(1, rowcount);
end
measurements.meas.m_b_unnormalized = mag_temp;
for i = 1: length(measurements.meas.m_b_unnormalized)
    if isnan(measurements.meas.m_b(1,i))
        measurements.meas.m_b_unnormalized(:,i) = measurements.meas.m_b_unnormalized(:,i-1);
        measurements.meas.m_b(:, i) = measurements.meas.m_b(:, i-1);
    end
end

measurements.available.gps = ones(1, rowcount);

i = find(contains(columns, 'timestampStartReadingSensors'));
timestampStartReadingSensors = M(i, :);
delta_t = [timestampStartReadingSensors, 0] - [0, timestampStartReadingSensors];
delta_t = delta_t(2:length(delta_t)-1);
if plotting
    plot(delta_t);
end

i = find(contains(columns, 'timestampFinishedReadingSensors'));
timestampFinishedReadingSensors = M(i, :);
delta_t = [timestampFinishedReadingSensors, 0] - [0, timestampFinishedReadingSensors];
delta_t = delta_t(2:length(delta_t)-1);

difference = timestampFinishedReadingSensors - timestampStartReadingSensors;

TA = median(delta_t)/1000; %[s]
TA = mean(delta_t)/1000; %[s]
time_s = linspace(0, TA*rowcount, rowcount);

measurements.meas.length_line = zeros(1, rowcount);
measurements.meas.gamma_line = zeros(2, rowcount);
           
end