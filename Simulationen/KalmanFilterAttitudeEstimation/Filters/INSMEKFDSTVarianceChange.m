classdef INSMEKFDSTVarianceChange < INSMEKF
    %% Multiplicative EKF INS with different sample time for the sensors
    % When the sensor is not available, the covariance is set to a high
    % value
    % Source: 
    % - [2]: Maley2013  - Multiplicative Quaternion Extended Kalman Filtering for Nonspinning Guided Projectiles
    % - [3]: Markley2003 - Attitude Error Representations for Kalman Filtering
    % - [4]: Patrick Spieler - https://github.com/Murmele/adcs
    % - [5]: Mattia Giurato - https://github.com/Murmele/Attitude_Estimation
    % AHRS in with quaternions
    % Input u: angular velocity in body coords, acceleration
    % Measurement y: magnetometer, position gps, height barometer, position
    % line angle sensor, velocity gps
    % y: [m_x, m_y, m_z, pxgps, pygps, pzgps, h_baro, pxline, pyline, pzline, vxgps, vygps, vzgps]';
    % x: - a: rotation error
    %    - p: position
    %    - v: velocity
    %    - delta b_w: error of bias
    %
    %
    % Difference to the INSMEKF:
    % - This filter is the same than the
    % "AHRSMEKFDifferentSampleTimeVarianceChange", but implemented as INS
    % and not as AHRS
    
    properties
        % are multiplied with sampleTime in initFilter function
        sigma_2_a_b_not_available = 0.1^2 * ones(3,1);
        sigma_2_yaw_not_available = 0.5^2 * ones(3,1);
        sigma_2_p_gps_not_available = 1 * ones(3,1);
        sigma_2_p_line_not_available = 1 * ones(3,1);
        sigma_2_v_gps_not_available = 1 * ones(3,1);
        sigma_2_h_baro_not_available = 1 * ones(3,1);
    end
    
    properties(Access = protected) % These variables must be initialised. Here or in the setupImpl function
        sigma2_process_noise;
        sigma2_measurement_noise;
        sigma2_measurement_noise_not_available;
        %
        v_gps_log;
        p_baro_log;
        p_gps_log;
        p_lineAngle_log;
    end
    
    methods
        
        function obj = INSMEKFDSTVarianceChange()
            obj@INSMEKF();
        end
        
        function initFilter(obj, sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, x0, P0, sigma2_process_noise, sigma2_measurement_noise)
            obj.sigma2_process_noise = sigma2_process_noise;
            obj.sigma2_measurement_noise = sigma2_measurement_noise; %/sampleTime;
            disp('INSMEKFDSTVarianceChange.m: Multiplying not available sensors with sigma_2*0.5');
            obj.sigma2_measurement_noise_not_available = sigma2_measurement_noise*0.5/sampleTime;
            initFilter@INSMEKF(obj, sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, x0, P0, diag(sigma2_process_noise), diag(sigma2_measurement_noise));
        end
        
        % estimate states
        function K = performEKF(obj, u, y, acc_available, gyro_available, mag_available, baro_available, gps_available, lineAngle_available)
            
            u(4:6) = u(4:6) * obj.g_mps2; % [g] --> [m/s^2]

            % interpolate obj.q from gyro
            adcs = 0;
            euler = 0;
            if ~adcs
                euler = 0;
                if euler
                    % derivation see EKFQuaternion.wxm
                    warning('Reimplement q_dot_f(u)!');
                    % do not work, because the indices are wrong. 
                    % reimplement this function!
                    q_dot = obj.q_dot_f(u); 
                    obj.q = obj.q + q_dot * obj.sampleTime; % estimate the new attitude
                    obj.q = obj.q / norm(obj.q); % normalize quaternion

                    b_w_dot = [0; 0; 0]; % always zero
                    obj.b_w = obj.b_w + b_w_dot * obj.sampleTime;
                else
                    gyro = u(1:3);
                    % [2], eq. 14
                    gyro_mean = (obj.gyro + gyro - obj.b_w)./2;

                    % [2], eq. 12
                    Omega = [0, -gyro_mean';
                             gyro_mean,      - obj.crossproductMatrix(gyro_mean)];
                    % [2], eq. 13     
                    obj.q = expm(1/2 * Omega * obj.sampleTime) * obj.q;
                    % needed to calculate integral of q_dot
                    obj.gyro = gyro - obj.b_w;
                end
            
            else
                % used from [4] 
                % MEKF3DGyro.m entnommen ( works same as the above approach)
                gyro = u(1:3);
                omega = gyro - obj.b_w; % omega_hat
                ang = norm(omega) * obj.sampleTime;
                if ang > 0.000001
                    axis = omega / norm(omega);
                    delta_q_ref = [cos(ang/2); axis*sin(ang/2)];  % where this eq. come from?
                else
                    delta_q_ref = [1; omega*obj.sampleTime/2];
                end
                q = Quaternion(obj.q);
                obj.q = q.mult(delta_q_ref);
            end
            
            [x_dot, G_lin] = obj.x_dot_and_G_lin(obj.EKF.x, u);
            obj.EKF.predictorStep(x_dot, G_lin);

            % wieso muss hier der reset gemacht werden?
            % Problem ist, dass dadurch, dass Beta nicht resetet wird, dies
            % delta a ungleich null setzt.
            % predicted state of the attitude error is zero.
            % Maybe the reset is here, because only alpha should be zero,
            % but in x_dot is not only an error state, when using all error
            % state there is no difference if the reset is here or after
            % the updateAttitude()
            obj.EKF.x = obj.EKF.x.*[0;0;0;1;1;1;1;1;1;1;1;1]; % reset [2]: eq. 92
                  
            available = [mag_available*ones(3,1); baro_available; gps_available*ones(3,1); lineAngle_available; gps_available*ones(3,1)];
            [y_hat, H_lin] = obj.y_hat_and_H_lin(obj.EKF.x, u, available);
            % integrate when not available
            if gps_available
                obj.p_gps_log = y(5:7);
                obj.v_gps_log = y(11:13);
            else
                p = obj.p_gps_log + obj.v_gps_log.*obj.sampleTime;
                y(5:7) = p;
            end
            if baro_available
                obj.p_baro_log = y(4);
            else
                p = obj.p_baro_log + obj.v_gps_log(3).*obj.sampleTime;
                y(4) = p;
            end
            if lineAngle_available
                obj.p_lineAngle_log = y(7:9);
            else
                p = obj.p_lineAngle_log + obj.v_gps_log.*obj.sampleTime;
                y(7:9) = p;
            end
            obj.EKF.correctorStep(y, y_hat, H_lin);
            
            obj.updateAttitude();
            % Wieso nicht hier den reset machen, sondern weiter oben?
            %obj.EKF.x = obj.EKF.x .* 0; % reset; this is wrong. works only, when all states are errors
            %obj.EKF.x = obj.EKF.x.*[0;0;0;1;1;1]; %reset only alpha
            
            K = obj.EKF.K;
        end
        
       % reimplementing xdot_and_G_lin for using also different sample
       % times

        % y_htt and H_lin were calculated in maxima:
        % "KalmanfilterLagebestimmungEKF_Systemmatrizen/Quaternionen.wxm"
        % Chapter: Standard Multiplicative Extended Kalman Filter
        
        function [y_hat, H_lin] = y_hat_and_H_lin(obj, x, u, available)
           delta_i = obj.delta_i;
            delta_d = obj.delta_d;
            g_mps2 = obj.g_mps2;

            % x: state vector 12x1
            a0 = x(1);
            a1 = x(2);
            a2 = x(3);
            px = x(4);
            py = x(5);
            pz = x(6);
            vx = x(7);
            vy = x(8);
            vz = x(9);
            b_w_phi = x(10); % bias
            b_w_theta = x(11);
            b_w_psi = x(12);

            w_phi_m= u(1);
            w_theta_m = u(2);
            w_psi_m = u(3);


            
            q0_hat = obj.q(1);
            q1_hat = obj.q(2);
            q2_hat = obj.q(3);
            q3_hat = obj.q(4);
            
            y_hat = [a2*(-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-a1*(-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat);
                    -a2*(cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat))-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))+a0*(-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat);
                    a1*(cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat))-a0*(-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat))-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2));
                    pz;
                    px;
                    py;
                    pz;
                    px;
                    py;
                    pz;
                    vx;
                    vy;
                    vz];
                    
           H_lin = [0,	2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)-2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)+sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)),	-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat),	0,	0,	0,	0,	0,	0,	0,	0,	0;
                    -2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)),	0,	-cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))+2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)+2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat),	0,	0,	0,	0,	0,	0,	0,	0,	0;
                    sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))+2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)-2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat),	cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0];
                
           % change variances, if the sensors are not available
           variances = obj.sigma2_measurement_noise;
           variances(~available) = obj.sigma2_measurement_noise_not_available(~available);
           obj.EKF.R_k = diag(variances);
        end
    end
end