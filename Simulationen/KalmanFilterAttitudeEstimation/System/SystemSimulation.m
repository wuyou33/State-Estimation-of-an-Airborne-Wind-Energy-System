% Simulation of the input values for the different filters
classdef SystemSimulation < handle
   properties
       sampleTime;
       TA_sensor; % sensor sample times
       delta_i; % inclination of magnetic field
       delta_d; % declination of magnetic field
       g_mps; % magnitude of gravity [m/s^2]
       sigma_w; % sigma_w^2: variance of gyro
       sigma_a; % variance of accelerometer
       sigma_b_w; % variance of gyro bias
       sigma_mag; % variance of magnetometer
       sigma_p_gps; % variance of position
       sigma_v_gps; % variance of velocity
       sigma_baro;
       sigma_p_line;
       sigma_gamma_line;
       sigma_length_line;

   end
   
   methods
        % a_i: acceleration in inertia coords
        % w_i: angular velocity in inertia coords
        function obj = SystemSimulation(sampleTime, TA_sensor, delta_i, delta_d, g_mps, sigma2)
            obj.sampleTime = sampleTime;
            obj.TA_sensor = TA_sensor;
            obj.g_mps = g_mps;
            obj.delta_i = delta_i;
            obj.delta_d = delta_d;
            
            obj.sigma_w = sqrt(sigma2.w_b);
            obj.sigma_a = sqrt(sigma2.a_b);
            obj.sigma_b_w = sqrt(sigma2.b_w);
            obj.sigma_mag = sqrt(sigma2.mag);
            obj.sigma_p_gps = sqrt(sigma2.p_gps);
            obj.sigma_v_gps = sqrt(sigma2.v_gps);
            obj.sigma_baro = sqrt(sigma2.baro);
            obj.sigma_p_line = sqrt(sigma2.p_line);
            obj.sigma_gamma_line = sqrt(sigma2.gamma_line);
            obj.sigma_length_line = sqrt(sigma2.length_line);
        end

        function [measurements, truthdata] = simulate(obj, timestamp, a_i, w_i, gamma_inertia0, p0, v0, b_w0, disableMagnetometer, sensorFailure)
             
           nrTime = length(timestamp);
            
           b_w_dot = [0;0;0];
           b_w = b_w0;
           v_i = zeros(3, nrTime);
           p_i = zeros(3, nrTime);
           v_i(:, 1) = v0;
           p_i(:, 1) = p0;
           gamma_inertia = gamma_inertia0;
           
           % calculate position and velocity
           i = 2:nrTime;
           v_i(:, i) = a_i(:, i) * obj.sampleTime; % euler integrator
           v_i = cumsum (v_i, 2);
           p_i(:, i) = v_i(:, i) * obj.sampleTime; % euler integrator
           p_i = cumsum (p_i, 2);
           
           a_i(3,:) = (a_i(3,:)-ones(1, nrTime).*obj.g_mps); % add gravity
           a_i = a_i/obj.g_mps;
           a_i(3,:) = a_i(3,:).* -1; % because z direction shows downward
           
           last_index_gps = 1;
           last_index_magnet = 1;
           last_index_acc = 1;
           last_index_gyro = 1;
           last_index_baro = 1;
           last_index_line = 1;
           p_meas = zeros(3, nrTime);
           v_meas = zeros(3, nrTime);
           a_b_meas = zeros(3, nrTime);
           w_b_meas = zeros(3, nrTime);
           m_b_meas = zeros(3, nrTime);
           baro_meas = zeros(1, nrTime);
           p_line_meas = zeros(3, nrTime);
           length_line_meas = zeros(1, nrTime);
           gamma_line_meas = zeros(2, nrTime);
           truthdata.gamma_inertia = zeros(3, nrTime);
           me = obj.Rmag_to_earth(obj.delta_i, obj.delta_d);
           
           gps_available = zeros(1, nrTime);
           acc_available = zeros(1, nrTime);
           gyr_available = zeros(1,nrTime);
           mag_available = zeros(1, nrTime);
           baro_available = zeros(1, nrTime);
           line_available = zeros(1, nrTime);
           for i=1: nrTime
               % gps (position and velocity)
               if (timestamp(i) - timestamp(last_index_gps) < obj.TA_sensor.gps && i > 1) || (sensorFailure.gps.fail == true && sensorFailure.gps.time <= timestamp(i))
                   if i == 1 % if sensor is from the beginning on not available
                       v_meas(:, i) = zeros(3,1);
                       p_meas(:, i) = zeros(3,1);
                   else
                        v_meas(:, i) = v_meas(:, i-1);
                        p_meas(:, i) = p_meas(:, i-1);
                   end
                   gps_available(i) = 0;
               else % v0, p0!
                   last_index_gps = i;
                   v_meas(:,i) = v_i(:, i) + obj.sigma_v_gps.*randn(3,1);
                   p_meas(:,i) = p_i(:, i) + obj.sigma_p_gps.*randn(3,1);
                   gps_available(i) = 1;
               end
               
               % barometer
               if (timestamp(i) - timestamp(last_index_baro) < obj.TA_sensor.baro && i > 1) || (sensorFailure.barometer.fail == true && sensorFailure.barometer.time <= timestamp(i))
                   if i == 1 % if sensor is from the beginning on not available
                        baro_meas(i) = 0;
                   else
                       baro_meas(i) = baro_meas(i-1);
                   end
                   baro_available(i) = 0;
               else
                   last_index_baro = i;
                   baro_meas(i) = p_i(3, i) + obj.sigma_baro*randn(1,1);
                   baro_available(i) = 1;
               end
               
               % line angle measurement
               if (timestamp(i) - timestamp(last_index_line) < obj.TA_sensor.line && i > 1) || (sensorFailure.lineAngle.fail == true && sensorFailure.lineAngle.time <= timestamp(i))
                   if i == 1 % if sensor is from the beginning on not available
                       p_line_meas(:, i) = p_line_meas(:, i-1);
                       gamma_line_meas(:, i) = gamma_line_meas(:, i-1);
                       length_line_meas(i) = length_line_meas(i-1);
                   else
                       p_line_meas(:, i) = zeros(3,1);
                       gamma_line_meas(:, i) = zeros(2,1);
                       length_line_meas(i) = 0;
                   end
                   line_available(i) = 0;
               else
                   last_index_line = i;
                   p_line_meas(:,i) = p_i(:, i) + obj.sigma_p_line.*randn(3,1);
                   theta = acos(p_i(3, i))/norm(p_i(:,i));
                   phi = atan2(p_i(2, i), p_i(1,i));
                   gamma_line_meas(:, i) = [phi; theta] + obj.sigma_gamma_line .* rand(2,1);
                   length_line_meas(i) = norm(p_i(:,i)) + obj.sigma_length_line .* rand(1,1);
                   line_available(i) = 1;
               end

               % acceleration in body coords
               if (timestamp(i) - timestamp(last_index_acc) < obj.TA_sensor.accelerometer && i > 1) || (sensorFailure.accelerometer.fail == true && sensorFailure.accelerometer.time <= timestamp(i))
                   if i == 1 % if sensor is from the beginning on not available
                       a_b_meas(:, i) = zeros(3,1);
                   else
                       a_b_meas(:, i) = a_b_meas(:, i-1);
                   end
                   acc_available(i) = 0;
               else
                   last_index_acc = i;
                   a_b = obj.rotate_from_earth_to_body_frame(a_i(:, i), gamma_inertia);
                   a_b_meas(:, i) = a_b + obj.sigma_a.*randn(3,1);
                   acc_available(i) = 1;
               end
               
               % angular velocity in body coords
               if (timestamp(i) - timestamp(last_index_gyro) < obj.TA_sensor.gyro && i > 1) || (sensorFailure.gyroscope.fail == true && sensorFailure.gyroscope.time <= timestamp(i))
                   if i == 1 % if sensor is from the beginning on not available 
                       w_b_meas(:, i) = zeros(3,1);
                   else
                       w_b_meas(:, i) = w_b_meas(:, i-1);
                   end
                    gyr_available(i) = 0;
               else
                   last_index_gyro = i;
                   w_b = obj.rotate_from_earth_to_body_frame(w_i(:, i), gamma_inertia);
                   % noise
                   w_b_meas(:, i) = w_b + obj.sigma_w .* randn(3,1) + b_w;
                   gyr_available(i) = 1;
               end
               
               % magnetometer
               if (timestamp(i) - timestamp(last_index_magnet) < obj.TA_sensor.magnetometer && i > 1) || (sensorFailure.magnetometer.fail == true && sensorFailure.magnetometer.time <= timestamp(i))
                   if i == 1 % if sensor is from the beginning on not available 
                       m_b_meas(:, i) = zeros(3,1);
                   else
                       m_b_meas(:, i) = m_b_meas(:, i-1);
                   end
                   mag_available(i) = 0;
               else
                   last_index_magnet = i;
                   m_b = obj.rotate_from_earth_to_body_frame(me, gamma_inertia);
                   m_b_meas(:, i) = m_b + obj.sigma_mag .* randn(3,1);
                   mag_available(i) = 1;
               end

               b_w = b_w + (b_w_dot+obj.sigma_b_w .* randn(3,1)) * obj.sampleTime;
               gamma_inertia = gamma_inertia + w_i(:, i) * obj.sampleTime;
               truthdata.gamma_inertia(:, i) = gamma_inertia;
           end
           
           if disableMagnetometer
                mag_available = zeros(1, nrTime);
           end
           
           measurements.meas.p_gps = p_meas;
           measurements.meas.v_gps = v_meas;
           measurements.meas.a_b = a_b_meas;
           measurements.meas.w_b = w_b_meas;
           measurements.meas.m_b = m_b_meas;
           measurements.meas.baro = baro_meas;
           measurements.meas.p_line = p_line_meas;
           measurements.meas.length_line = length_line_meas;
           measurements.meas.gamma_line = gamma_line_meas;
           measurements.available.gps = gps_available;
           measurements.available.acc = acc_available;
           measurements.available.mag = mag_available;
           measurements.available.gyr = gyr_available;
           measurements.available.baro = baro_available;
           measurements.available.line = line_available;
           truthdata.v_i = v_i;
           truthdata.p_i = p_i;
                 
           
        end
   end
   
   methods (Static)
       function y = rotate_from_earth_to_body_frame(v,gamma)
            % variable re_to_b in Systemmatrix.wxmx
            phi = gamma(1);
            theta = gamma(2);
            psi = gamma(3);
            y= 	[cos(psi)*cos(theta),	sin(psi)*cos(theta),	-sin(theta);
                sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi),	sin(phi)*sin(psi)*sin(theta)+cos(phi)*cos(psi),	sin(phi)*cos(theta);
                cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi),	cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi),	cos(phi)*cos(theta)]*v;
       end
        
       function me = magentic_field(delta_i, delta_d)
           % from mag to earth
           me=     [cos((pi*delta_d)/180)*cos((pi*delta_i)/180);
                     -sin((pi*delta_d)/180)*cos((pi*delta_i)/180);
                     -sin((pi*delta_i)/180)];
       end
       
       function R = Re_to_mag(delta_i, delta_d)
           R = [cos((pi*delta_d)/180)*cos((pi*delta_i)/180),	-sin((pi*delta_d)/180)*cos((pi*delta_i)/180),	-sin((pi*delta_i)/180);
                sin((pi*delta_d)/180),	cos((pi*delta_d)/180),	0;
                cos((pi*delta_d)/180)*sin((pi*delta_i)/180),	-sin((pi*delta_d)/180)*sin((pi*delta_i)/180),	cos((pi*delta_i)/180)];
       end
       
       function v = Rmag_to_earth(delta_i, delta_d)
           v = [cos((pi*delta_d)/180)*cos((pi*delta_i)/180);
                -sin((pi*delta_d)/180)*cos((pi*delta_i)/180);
                -sin((pi*delta_i)/180)];
       end
   end
   
end