classdef AHRSMEKFDifferentSampleTimeMagnetometer < AHRSMEKF
    %% Multiplicative EKF AHRS
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
    % - Magnetometer is used only for estimating the yaw angle and not the
    % magnetic vector is used completely in the kalman filter
    properties
    end
    
    methods
        
        function obj = AHRSMEKFDifferentSampleTimeMagnetometer(sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw)
            obj@AHRSMEKF();
            obj.initFilter(sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw)
        end
        
        function initFilter(obj, sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw)
            initFilter@AHRSMEKF(obj, sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw);
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
                  
            available = [acc_available, mag_available];
            for i=0:1
                if available(i+1)
                    y_part = y;%(i*3 +1 : i*3+1+2);
                    [y_meas, y_hat, H_lin] = obj.y_hat_and_H_lin(obj.EKF.x, u, y_part, i);
                    obj.EKF.correctorStep(y_meas, y_hat, H_lin);
                end
            end
            obj.updateAttitude();
            % Wieso nicht hier den reset machen, sondern weiter oben?
            %obj.EKF.x = obj.EKF.x .* 0; % reset; this is wrong. works only, when all states are errors
            %obj.EKF.x = obj.EKF.x.*[0;0;0;1;1;1]; %reset only alpha
            
            K = obj.EKF.K;
        end

        % y_htt and H_lin were calculated in maxima:
        % "KalmanfilterLagebestimmungEKF_Systemmatrizen/Quaternionen.wxm"
        % Chapter: Standard Multiplicative Extended Kalman Filter; Measurement matrix (Use magnetometer only for yaw angle estimation)
        function [y_meas, y_hat, H_lin] = y_hat_and_H_lin(obj, x, u, y, nr)
            
            if nr == 0
                % accelerometer
                q0_hat = obj.q(1);
                q1_hat = obj.q(2);
                q2_hat = obj.q(3);
                q3_hat = obj.q(4);
                g_mps2 = obj.g_mps2;
                a0 = x(1);
                a1 = x(2);
                a2 = x(3);
                b_w_phi = x(4); % bias
                b_w_theta = x(5);
                b_w_psi = x(6);
                v_x = sqrt(u(4)^2+u(5)^2);
                v_y = 0;
                v_z = 0;
                w_phi_m= u(1);
                w_theta_m = u(2);
                w_psi_m = u(3);
                
                y_meas = y(1:3);

                H_lin =  [0,	-g_mps2*(1-2*(q2_hat^2+q1_hat^2)),	2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat),	0,	-v_z,	v_y;
                            g_mps2*(1-2*(q2_hat^2+q1_hat^2)),	0,	-2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat),	v_z,	0,	-v_x;
                            -2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat),	2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat),	0,	-v_y,	v_x,	0];

                y_hat = [v_z*(w_theta_m-b_w_theta)+v_y*(b_w_psi-w_psi_m)+2*a2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat)+2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat)-a1*g_mps2*(1-2*(q2_hat^2+q1_hat^2));
                        v_x*(w_psi_m-b_w_psi)+v_z*(b_w_phi-w_phi_m)+2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat)-2*a2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat)+a0*g_mps2*(1-2*(q2_hat^2+q1_hat^2));
                        v_x*(b_w_theta-w_theta_m)+v_y*(w_phi_m-b_w_phi)-2*a0*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat)+2*a1*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat)+g_mps2*(1-2*(q2_hat^2+q1_hat^2))];
                obj.EKF.R_k = obj.R(1:3,1:3)/ obj.sampleTime;
            else
                % magnetometer
                q0_hat = obj.q(1);
                q1_hat = obj.q(2);
                q2_hat = obj.q(3);
                q3_hat = obj.q(4);
                a0 = x(1);
                a1 = x(2);
                a2 = x(3);
                g_mps2 = obj.g_mps2;
                v_x = sqrt(u(4)^2+u(5)^2);
                v_y = 0;
                v_z = 0;
                w_phi_m= u(1);
                w_theta_m = u(2);
                w_psi_m = u(3);
                gyro = [w_phi_m;
                        w_theta_m;
                        w_psi_m];
                mag = y(4:6);
                
                q_temp = Quaternion(obj.q);
                alpha = [a0; a1; a2];
                qAlpha = [2; alpha];
                q_hat2 = q_temp.mult(qAlpha);
                q_hat2= q_hat2/norm(q_hat2);
                q0_hat2 = q_hat2(1);
                q1_hat2 = q_hat2(2);
                q2_hat2 = q_hat2(3);
                q3_hat2 = q_hat2(3);
                
                v_vel_b = [sqrt(v_x^2+v_y^2); 0;0];
                % must be multiplied by -1 I think, because the accelerometer
                % measures the centrifugal acceleration??
                v_acc_centripetal = obj.crossproductMatrix(gyro - obj.EKF.x(4:6))*v_vel_b; % ignored in the original paper.
                
                % it is unstable to calculate acc_hat from estimated q,
                % with delta_d and delta_i > 0 and gamma_inertia0 = [0; 0; pi/4];
                % with lower acceleration noise it is better
                acc_hat = [2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat);
                            2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat);
                            g_mps2*(1-2*(q2_hat^2+q1_hat^2))];
                acc_hat2 = [2*g_mps2*(q1_hat2*q3_hat2-q0_hat2*q2_hat2);
                            2*g_mps2*(q2_hat2*q3_hat2+q0_hat2*q1_hat2);
                            g_mps2*(1-2*(q2_hat2^2+q1_hat2^2))];
                acc_hat_meas = y(1:3); %  - v_acc_centripetal; % centipetal compensation
                [phi_meas3, theta_meas3, y_meas3] = calculateEulerAnglesFromAccelerometerAndMagnetometer(acc_hat, mag, obj.delta_d, obj.delta_i);
                [phi_meas2, theta_meas2, y_meas2] = calculateEulerAnglesFromAccelerometerAndMagnetometer(acc_hat2, mag, obj.delta_d, obj.delta_i);
                % only psi is interesting here
                [phi_meas, theta_meas, y_meas] = calculateEulerAnglesFromAccelerometerAndMagnetometer(acc_hat_meas, mag, obj.delta_d, obj.delta_i);
                y_meas = deg2rad(y_meas);

                H_lin = [(2*(-(q2_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+(q3_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+(q0_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-(q1_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)*(-2*q1_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)+2*q0_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)+2*q3_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)-2*q2_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2))-((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)*(-2*q1_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)+2*q0_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)+2*q3_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)-2*q2_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2)))*(1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))))/((1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))^(2)+4*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))^(2))+(4*((2*q3_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-(2*q2_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2*(-2*q1_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)+2*q0_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)+2*q3_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)-2*q2_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2))-((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2*(-2*q1_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)+2*q0_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)+2*q3_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)-2*q2_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2)))*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))/((1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))^(2)+4*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))^(2)),(2*((q1_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+(q0_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-(q3_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-(q2_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)*(-2*q2_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)-2*q3_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)+2*q0_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)+2*q1_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2))-((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)*(-2*q2_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)-2*q3_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)+2*q0_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)+2*q1_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2)))*(1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))))/((1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))^(2)+4*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))^(2))+(4*((2*q0_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+(2*q1_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2*(-2*q2_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)-2*q3_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)+2*q0_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)+2*q1_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2))-((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2*(-2*q2_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)-2*q3_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)+2*q0_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)+2*q1_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2)))*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))/((1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))^(2)+4*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))^(2)),(2*((q0_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-(q1_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+(q2_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-(q3_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)*(-2*q3_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)+2*q2_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)-2*q1_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)+2*q0_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2))-((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)*(-2*q3_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)+2*q2_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)-2*q1_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)+2*q0_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2)))*(1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))))/((1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))^(2)+4*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))^(2))+(4*(-(2*q1_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+(2*q0_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)-((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2*(-2*q3_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)+2*q2_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)-2*q1_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)+2*q0_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2))-((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2*(-2*q3_hat*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)+2*q2_hat*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)-2*q1_hat*(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)+2*q0_hat*(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)))/(((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)^(2)))*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))/((1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)))^(2)+4*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2))^(2)),0,0,0];
                y_hat = atan2(2*(((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)*(-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)*(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat))/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)),1-2*(((a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)+((2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)/((-a2*q3_hat-a1*q2_hat-a0*q1_hat+2*q0_hat)^2+(-a1*q3_hat+a2*q2_hat+2*q1_hat+a0*q0_hat)^2+(a0*q3_hat+2*q2_hat-a2*q1_hat+a1*q0_hat)^2+(2*q3_hat-a0*q2_hat+a1*q1_hat+a2*q0_hat)^2)));
                obj.EKF.R_k = obj.R(4, 4) / obj.sampleTime; % because only yaw is used
                
                % Did not solve the problem!
                % Problem: if one angle becomes greater than 180° ther is a
                % wrap around to -180° (pi rad). In this short time the angle
                % difference is 360°!              
                diff = pi - abs(abs(y_meas - y_hat) - pi);
                if y_meas < y_hat
                    diff = diff *-1;
                end
%                 y_meas = 0;
%                 y_hat = -diff;
                
            end
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