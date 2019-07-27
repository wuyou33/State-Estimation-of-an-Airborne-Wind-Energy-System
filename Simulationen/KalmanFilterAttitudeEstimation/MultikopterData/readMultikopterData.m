function [measurements, truthdata] = readMultikopterData(file)

if strcmp(file, '')
    error('No file specified!');
end

gyx = ncread(file, 'input_bus.mpu6050data.gyx')';
gyy = ncread(file, 'input_bus.mpu6050data.gyy')';
gyz = ncread(file, 'input_bus.mpu6050data.gyz')';
mpu6050_available = ncread(file, 'input_bus.mpu6050data.mpu6050_available')';
measurements.meas.w_b = [gyx;
                 gyy;
                 gyz] * pi / 180;
measurements.available.w_b = mpu6050_available;
cx = ncread(file, 'multikopter_signals.cx_model')';
cy = ncread(file, 'multikopter_signals.cy_model')';
cz = ncread(file, 'multikopter_signals.cz_model')';
hmc5883_available = ncread(file, 'input_bus.hmc5883data.hmc5883_available')';
measurements.meas.m_b = [cx;
                cy;
                cz];
measurements.meas.m_b_unnormalized = measurements.meas.m_b;
measurements.meas.m_b = measurements.meas.m_b ./ vecnorm(measurements.meas.m_b);
measurements.available.mag = hmc5883_available;

ax = ncread(file, 'multikopter_signals.acx')';
ay = ncread(file, 'multikopter_signals.acy')';
az = ncread(file, 'multikopter_signals.acz')';
measurements.meas.a_b = [ax;
                 ay;
                 az];
measurements.available.acc = mpu6050_available;
roll_deg = ncread(file, 'multikopter_signals.roll_deg')';
pitch_deg = ncread(file, 'multikopter_signals.pitch_deg')';
yaw_deg = ncread(file, 'multikopter_signals.yaw_deg')';
truthdata.gamma_inertia = [roll_deg;
                    pitch_deg;
                    yaw_deg];
truthdata.multikopterEstimatioGammaInertia = [roll_deg; pitch_deg; yaw_deg];
                
%% GPS                
lat = ncread(file, 'input_bus.gpsdata.latitude')';
long = ncread(file, 'input_bus.gpsdata.longitude')';
alt = ncread(file, 'input_bus.gpsdata.altitude')';
measurements.available.gps = mpu6050_available;
measurements.meas.p_gps = [lat; long; alt];
measurements.meas.v_gps = zeros(3, length(mpu6050_available));

%% Barometer
measurements.meas.baro = ncread(file, 'input_bus.bmp180data.altitude')';
measurements.available.baro = mpu6050_available;

time = ncread(file, 'time')';
TA = time(2)-time(1); % [µs]
TA = TA*1e-6;
f = 1/TA;
time = 0:TA:TA*(length(time)-1);
measurements.meas.time = time;
measurements.meas.TA = TA; % [s]

end
