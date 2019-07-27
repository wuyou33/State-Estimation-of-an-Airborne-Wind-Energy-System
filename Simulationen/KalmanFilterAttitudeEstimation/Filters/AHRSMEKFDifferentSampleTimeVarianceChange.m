classdef AHRSMEKFDifferentSampleTimeVarianceChange < AHRSMEKF
    %% Multiplicative EKF AHRS with different sample time for the sensors
    % Source: 
    % - [2]: Maley2013  - Multiplicative Quaternion Extended Kalman Filtering for Nonspinning Guided Projectiles
    % - [3]: Markley2003 - Attitude Error Representations for Kalman Filtering
    % - [4]: Patrick Spieler - https://github.com/Murmele/adcs
    % - [5]: Mattia Giurato - https://github.com/Murmele/Attitude_Estimation
    % AHRS in with quaternions
    % Input u: angular velocity in body coords
    % Measurement y: acceleration in body coords, magnetometer measurements
    % y: [acc_x, acc_y, acc_z, m_x, m_y, m_z]';
    % x: - a: rotation error
    %    - delta b_w: error of bias
    %
    %
    % Difference to the AHRSMEKFDifferentSampleTime:
    % - the correction step is calculated only once, but if the measurement
    % is not available, the measurement matrix entries are changed to
    % higher values to simulate that the sensor is not available
    
    % Problems:
    % - When magnetometer is always off, the filter is unstable!
    % -- Solution: Use same principle like in
    % "AHRSMEKFDifferentSampleTimeMagnetometer2" to set the
    % difference of the other sensor to zero and use the original equations
    properties
        % are multiplied with sampleTime in initFilter function
        sigma_2_a_b_not_available = 0.1^2 * ones(3,1);
        sigma_2_yaw_not_available = 0.5^2 * ones(3,1);
    end
    
    methods
        
        function obj = AHRSMEKFDifferentSampleTimeVarianceChange(sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw)
            obj@AHRSMEKF();
            obj.initFilter(sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw)
        end
        
        function initFilter(obj, sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw)
            initFilter@AHRSMEKF(obj, sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw);
            obj.sigma_2_a_b_not_available = obj.sigma_2_a_b_not_available * sampleTime;
            obj.sigma_2_yaw_not_available = obj.sigma_2_yaw_not_available * sampleTime;
        end
        
        % estimate states
        function K = performEKF(obj, u, y, acc_available, mag_available)
            
            y(1:3) = y(1:3) * obj.g_mps2; % [g] --> [m/s^2]

            % interpolate obj.q from gyro
            adcs = 0;
            euler = 0;
            if ~adcs
                euler = 0;
                if euler
                    % derivation see EKFQuaternion.wxm
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
                omega = gyro - obj.EKF.x(4:6); % omega_hat
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
            obj.EKF.x = obj.EKF.x.*[0;0;0;1;1;1]; % reset [2]: eq. 92
                  
            available = [acc_available*ones(3,1); mag_available*ones(3,1)];
            [y_hat, H_lin] = obj.y_hat_and_H_lin(obj.EKF.x, u, available);
            obj.EKF.correctorStep(y, y_hat, H_lin);
            
            obj.updateAttitude();
            % Wieso nicht hier den reset machen, sondern weiter oben?
            %obj.EKF.x = obj.EKF.x .* 0; % reset; this is wrong. works only, when all states are errors
            %obj.EKF.x = obj.EKF.x.*[0;0;0;1;1;1]; %reset only alpha
            
            K = obj.EKF.K;
        end

        % y_htt and H_lin were calculated in maxima:
        % "KalmanfilterLagebestimmungEKF_Systemmatrizen/Quaternionen.wxm"
        % Chapter: Standard Multiplicative Extended Kalman Filter
        function [y_hat, H_lin] = y_hat_and_H_lin(obj, x, u, available)
                        delta_i = obj.delta_i;
            delta_d = obj.delta_d;
            g_mps2 = obj.g_mps2;

            % x: state vector 6x1
            a0 = x(1);
            a1 = x(2);
            a2 = x(3);
            b_w_phi = x(4); % bias
            b_w_theta = x(5);
            b_w_psi = x(6);

            w_phi_m= u(1);
            w_theta_m = u(2);
            w_psi_m = u(3);

            % in earth frame
            vx_ins = u(4);
            vy_ins = u(5);
            vz_ins = u(6);
            
            v_x = sqrt(u(4)^2+u(5)^2);
            v_y = 0;
            v_z = 0;
            
            q0_hat = obj.q(1);
            q1_hat = obj.q(2);
            q2_hat = obj.q(3);
            q3_hat = obj.q(4);

            H_lin = [0,	-g_mps2*(1-2*(q2_hat^2+q1_hat^2)),	2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat),	0,	-v_z,	v_y;
                    g_mps2*(1-2*(q2_hat^2+q1_hat^2)),	0,	-2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat),	v_z,	0,	-v_x;
                    -2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat),	2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat),	0,	-v_y,	v_x,	0;
                    0,	2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)-2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)+sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)),	-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat),	0,	0,	0;
                    -2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)),	0,	-cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))+2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)+2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat),	0,	0,	0;
                    sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))+2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)-2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat),	cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat),	0,	0,	0,	0];

           % identity matrix included
           y_hat = [v_z*(w_theta_m-b_w_theta)+v_y*(b_w_psi-w_psi_m)+2*a2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat)+2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat)-a1*g_mps2*(1-2*(q2_hat^2+q1_hat^2));
                    v_x*(w_psi_m-b_w_psi)+v_z*(b_w_phi-w_phi_m)+2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat)-2*a2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat)+a0*g_mps2*(1-2*(q2_hat^2+q1_hat^2));
                    v_x*(b_w_theta-w_theta_m)+v_y*(w_phi_m-b_w_phi)-2*a0*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat)+2*a1*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat)+g_mps2*(1-2*(q2_hat^2+q1_hat^2));
                    a2*(-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat))+cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-a1*(-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat);
                    -a2*(cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat))-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))+a0*(-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat);
                    a1*(cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat))-a0*(-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat))-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2))];

        
           % change variances, if the sensors are not available
           variances = [obj.sigma_2_a; obj.sigma_2_yaw];
           variances(~available) = [obj.sigma_2_a_b_not_available(~available(1:3)); obj.sigma_2_yaw_not_available(~available(4:6))];
           obj.EKF.R_k = diag(variances');
        end
    end

    methods (Static)
        function M = crossproductMatrix(v)
            M = [0,-v(3),v(2); 
                 v(3),0,-v(1);
                 -v(2),v(1),0];
        end
    end

end