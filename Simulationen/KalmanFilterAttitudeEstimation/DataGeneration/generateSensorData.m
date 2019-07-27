% source: Markus Sedlmair
% Markus did not use different sample rates
function [measurements, truthdata] = generateSensorData(kiteTruth, TA, TA_sensor, g_mps2, delta_i, delta_d, p0, v0, b_w0, sigma2, disableGPSatToHighAccelerations, tetherSag, sensorFailure)

% TA_sensor.gyro not used, because in the kiteTruth structure is already the
% angular rate

nbrOfDatas = length(kiteTruth.time_s);

measurements = [];

% attitude
yaw_rad   = pi/180*kiteTruth.eul_deg(3,:);
pitch_rad = pi/180*kiteTruth.eul_deg(2,:);
roll_rad  = pi/180*kiteTruth.eul_deg(1,:);
gamma_inertia = [roll_rad; 
                 pitch_rad;
                 yaw_rad];

% bias of gyro
b_w_dot = [0;0;0];
b_n = sigma2.b_w.*randn(3,nbrOfDatas); % noise of bias of gyro
b_w = zeros(3, nbrOfDatas);
b_w(:,1) = b_w0;
for i=2:nbrOfDatas
    b_w(:, i) = b_w(:, i-1) + b_w_dot*TA;
end
b_w = b_w + b_n;

fn = fieldnames(sensorFailure);
minTimeSensorFailure = 0;
sensorFails = 0;
first = 1;
for k=1:numel(fn)
    if strcmp(fn{k}, 'landing')
        continue;
    end
    if (sensorFailure.(fn{k}).time < minTimeSensorFailure || first == 1) && sensorFailure.(fn{k}).fail == true
        minTimeSensorFailure = sensorFailure.(fn{k}).time;
        sensorFails = 1;
        first = 0;
    end
end
minIndexSensorFailure = find(kiteTruth.time_s >= minTimeSensorFailure & logical(sensorFails*ones(1, nbrOfDatas)));
if ~isempty(minIndexSensorFailure)
    kiteTruth.pe_m(:, minIndexSensorFailure) = kiteTruth.pe_m(:, minIndexSensorFailure(1))*ones(1, length(minIndexSensorFailure));
    kiteTruth.ve_mps(:, minIndexSensorFailure) = zeros(3, length(minIndexSensorFailure));
    kiteTruth.ae_mps2(:, minIndexSensorFailure) = zeros(3, length(minIndexSensorFailure));
    kiteTruth.wb_rps(:, minIndexSensorFailure) = zeros(3, length(minIndexSensorFailure));
    kiteTruth.eul_deg(:, minIndexSensorFailure) = zeros(3, length(minIndexSensorFailure));
    gamma_inertia(:, minIndexSensorFailure) = zeros(3, length(minIndexSensorFailure));
end
clear fn minTimeSensorFailure first;
%% angular velocity in body coords
measurements.meas.w_b = zeros(3, nbrOfDatas);
measurements.available.gyr = zeros(1, nbrOfDatas);
last_index = 1;
for i = 1 : nbrOfDatas
    if (kiteTruth.time_s(i) - kiteTruth.time_s(last_index) < TA_sensor.gyro && i > 1) || (sensorFailure.gyroscope.fail == true && sensorFailure.gyroscope.time <= kiteTruth.time_s(i))
        if i == 1 % if sensor is from the beginning on not available
            measurements.meas.w_b(:,i) = zeros(3,1);
        else
            measurements.meas.w_b(:,i) = measurements.meas.w_b(:,i-1);
        end
        measurements.available.gyr(i) = 0;
        continue;
    end
    last_index = i;
    measurements.available.gyr(i) = 1;
    measurements.meas.w_b(:, i) = kiteTruth.wb_rps(:, i) + b_w(:, i) + sqrt(sigma2.w_b) .* randn(3, 1);
end

%% acceleration in body coords
gravity = [0; 0; g_mps2];
a_i = (kiteTruth.ae_mps2+gravity);
%a_i(3, :) = a_i(3, :) * -1; % because z direction shows downward
measurements.available.acc = zeros(1, nbrOfDatas);
last_index = 1;
for i=1: nbrOfDatas
    if (kiteTruth.time_s(i) - kiteTruth.time_s(last_index) < TA_sensor.accelerometer && i > 1) || (sensorFailure.accelerometer.fail == true && sensorFailure.accelerometer.time <= kiteTruth.time_s(i))
        if i == 1 % if sensor is from the beginning on not available
            measurements.meas.a_b(:,i) = zeros(3,1);
        else
            measurements.meas.a_b(:,i) = measurements.meas.a_b(:,i-1);
        end
        measurements.available.acc(i) = 0;
        continue;
    end
    measurements.available.acc(i) = 1;
    last_index = i;
    a_b = SystemSimulation.rotate_from_earth_to_body_frame(a_i(:, i), gamma_inertia(:, i)); % centrifugal forces already included
    a_b = a_b/g_mps2; % [m/s^2] --> [g]
    measurements.meas.a_b(:,i) = a_b + sqrt(sigma2.a_b).*randn(3, 1);
end

measurements.available.gps(vecnorm(measurements.meas.a_b)>3.5) = 0;

%% position and velocity
measurements.available.gps = zeros(1, nbrOfDatas);
measurements.meas.p_gps = zeros(3, nbrOfDatas);
measurements.meas.v_gps = zeros(3, nbrOfDatas);
last_index = 1;
gpsThresholdAcc = 4;
if disableGPSatToHighAccelerations
    disp(['The threshold maximum acceleration, the GPS module works, is set to ' num2str(gpsThresholdAcc) ' g']);
end
vecnorm_a_b = vecnorm(measurements.meas.a_b);
for i=1: nbrOfDatas
    
    if (kiteTruth.time_s(i) - kiteTruth.time_s(last_index) < TA_sensor.gps && i > 1) || (vecnorm_a_b(i) >gpsThresholdAcc && disableGPSatToHighAccelerations) || (sensorFailure.gps.fail == true && sensorFailure.gps.time <= kiteTruth.time_s(i))
        if i == 1 % if sensor is from the beginning on not available
            measurements.meas.v_gps(:, i) = zeros(3,1);
            measurements.meas.p_gps(:, i) = zeros(3,1);
        else
            measurements.meas.v_gps(:, i) = measurements.meas.v_gps(:, i-1);
            measurements.meas.p_gps(:, i) = measurements.meas.p_gps(:, i-1);
        end
        measurements.available.gps(i) = 0;
        continue;
    end
    measurements.available.gps(i) = 1;
    last_index = i;
    if i > 1
        measurements.meas.v_gps(:, i) = kiteTruth.ve_mps(:, i) + sqrt(sigma2.v_gps).*randn(3,1);
        x = measurements.meas.p_gps(:, i-1) + measurements.meas.v_gps(:, i) * TA_sensor.gps; % simple euler integrator
        measurements.meas.p_gps(:, i) = kiteTruth.pe_m(:, i) + sqrt(sigma2.p_gps).*randn(3,1);
    else
        measurements.meas.v_gps(:, i) = kiteTruth.ve_mps(:, i) + sqrt(sigma2.v_gps).*randn(3,1);
        measurements.meas.p_gps(:, i) = kiteTruth.pe_m(:, i) + sqrt(sigma2.p_gps).*randn(3,1);
    end
end

%% barometer
measurements.available.baro = zeros(1, nbrOfDatas);
measurements.meas.baro = zeros(1, nbrOfDatas);
last_index = 1;
for i=1: nbrOfDatas
    if (kiteTruth.time_s(i) - kiteTruth.time_s(last_index) < TA_sensor.baro && i > 1)  || (sensorFailure.barometer.fail == true && sensorFailure.barometer.time <= kiteTruth.time_s(i))
        if i == 1 % if sensor is from the beginning on not available
            measurements.meas.baro(i) = zeros(3,1);
        else
            measurements.meas.baro(i) = measurements.meas.baro( i-1);
        end
        measurements.available.baro(i) = 0;
        continue;
    end
    measurements.available.baro(i) = 1;
    last_index = i;
    measurements.meas.baro(i) = kiteTruth.pe_m(3, i) + sqrt(sigma2.baro).*randn(1);
end

%% line angle
v = kiteTruth.ve_mps(:, 1);
x = kiteTruth.pe_m(:, 1);
measurements.available.line = zeros(1, nbrOfDatas);
measurements.meas.p_line = zeros(3, nbrOfDatas);
% assumption: straight line
measurements.meas.gamma_line = zeros(2, nbrOfDatas);
measurements.meas.length_line = zeros(1, nbrOfDatas);

% assumption: line has sag due to the gravity
measurements.meas.theta_base_line = zeros(1, nbrOfDatas);
measurements.meas.theta_kite_line = zeros(1, nbrOfDatas);
measurements.meas.phi_sag_line = zeros(1, nbrOfDatas);
measurements.meas.length_sag_line = zeros(1, nbrOfDatas); % real length of the tether
last_index = 1;

l = 155;
disp(['Assumption: Tetherlength is always ', num2str(l), 'm']);
lerror = l + 0; % [m] measured length with error
disp(['Assumption: Measured tetherlenght has an error of ', num2str(lerror- l), 'm']);

for i=1: nbrOfDatas
    if (kiteTruth.time_s(i) - kiteTruth.time_s(last_index) < TA_sensor.line && i > 1) || (sensorFailure.lineAngle.fail == true && sensorFailure.lineAngle.time <= kiteTruth.time_s(i))
        if i == 1 % if sensor is from the beginning on not available
            measurements.meas.p_line(:, i) = zeros(3,1);
            measurements.meas.gamma_line(:, i) = zeros(2,1);
            measurements.meas.length_line(i) = 0;

            measurements.meas.theta_base_line(i) = 0;
            measurements.meas.theta_kite_line(i) = 0;
            measurements.meas.phi_sag_line(i) = 0;
            measurements.meas.length_sag_line(i) = 0;
        else
            measurements.meas.p_line(:, i) = measurements.meas.p_line(:, i-1);
            measurements.meas.gamma_line(:, i) = measurements.meas.gamma_line(:, i-1);
            measurements.meas.length_line(i) = measurements.meas.length_line(i -1);

            measurements.meas.theta_base_line(i) = measurements.meas.theta_base_line(i-1);
            measurements.meas.theta_kite_line(i) = measurements.meas.theta_kite_line(i-1);
            measurements.meas.phi_sag_line(i) = measurements.meas.phi_sag_line(i -1);
            measurements.meas.length_sag_line(i) = measurements.meas.length_sag_line(i-1);
        end
        measurements.available.line(i) = 0;
        continue;
    end
    measurements.available.line(i) = 1;
    last_index = i;
    
    if ~tetherSag
        l = norm(kiteTruth.pe_m(:,i));
    end
    theta = acos(max(kiteTruth.pe_m(3, i),0)/l); % don't allow values smaller than zero
    phi = atan2(kiteTruth.pe_m(2, i) , kiteTruth.pe_m(1,i));
    measurements.meas.gamma_line(:, i) = [theta; phi] + sqrt(sigma2.gamma_line).*randn(2,1);
    
    measurements.meas.length_line(i) = l + sqrt(sigma2.length_line).*randn(1);
    
    if tetherSag
        % sag in tether
        pos = kiteTruth.pe_m(:, i) + sqrt(sigma2.p_line).*randn(3,1);
        [theta_base, theta_kite, phi_sag, ~, ~] = calculateAnglesFromPosition(pos(1), pos(2), pos(3),l);
        measurements.meas.theta_base_line(i) = theta_base + sqrt(sigma2.gamma_line(1))*randn(1);
        measurements.meas.theta_kite_line(i) = theta_kite + sqrt(sigma2.gamma_line(1))*randn(1);
        measurements.meas.phi_sag_line(i) = phi_sag + sqrt(sigma2.gamma_line(1))*randn(1);
        measurements.meas.length_sag_line(i) = lerror + sqrt(sigma2.length_line) * randn(1); % assumption: 
    else
        measurements.meas.theta_base_line(i) =  measurements.meas.gamma_line(1, i);
        measurements.meas.theta_kite_line(i) =  measurements.meas.gamma_line(1, i);
        measurements.meas.phi_sag_line(i) = measurements.meas.gamma_line(2, i);
        measurements.meas.length_sag_line(i) = measurements.meas.length_line(i);
    end
    
    if tetherSag
        disp('"measurements.meas.p_line" is manipulated here');
        warning('With this manipulation some filters do not work pretty well!');
        p_line_calc = [ lerror.*sin(measurements.meas.theta_base_line(i)) .* cos(measurements.meas.phi_sag_line(i)); 
                        lerror.*sin(measurements.meas.theta_base_line(i)) .* sin(measurements.meas.phi_sag_line(i));
                        lerror.*cos(measurements.meas.theta_base_line(i))];
        measurements.meas.p_line(:, i) = p_line_calc; % + sqrt(sigma2.p_line).*randn(3,1);
    else
        measurements.meas.p_line(:, i) = kiteTruth.pe_m(:, i) + sqrt(sigma2.p_line).*randn(3,1);
    end
end
clear theta_base theta_kite phi_sag

%% magnetometer
me=     Rmag_to_e([1;0;0], delta_d, delta_i);
last_index = 1;
measurements.available.mag = zeros(1, nbrOfDatas);
for i = 1: nbrOfDatas
    if (kiteTruth.time_s(i) - kiteTruth.time_s(last_index) < TA_sensor.magnetometer && i > 1)  || (sensorFailure.magnetometer.fail == true && sensorFailure.magnetometer.time <= kiteTruth.time_s(i))
        if i == 1 % if sensor is from the beginning on not available
            measurements.meas.m_b(:, i) = zeros(3,1);
        else
            measurements.meas.m_b(:, i) = measurements.meas.m_b(:, i-1);
        end
        measurements.available.mag(i) = 0;
        continue;
    end
    measurements.available.mag(i) = 1;
    last_index = i;
    m_b = SystemSimulation.rotate_from_earth_to_body_frame(me, gamma_inertia(:, i));
    measurements.meas.m_b(:, i) = m_b + sqrt(0).*randn(3, 1); % TODO: set noise value
end

truthdata.p_i = kiteTruth.pe_m;
truthdata.v_i = kiteTruth.ve_mps;
truthdata.gamma_inertia = kiteTruth.eul_deg;

end