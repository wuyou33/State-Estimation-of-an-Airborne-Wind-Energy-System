classdef AHRSMEKF < matlab.System
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
    % Version: 1.0
    
    % These variables must be initialized in the simulink model
    properties
        sigma_2_w;
        sigma_2_a;
        sigma_2_b_w;
        sigma_2_yaw;
        sampleTime;
        delta_i;
        delta_d;
        g_mps2;
    end
    
    properties(Access = protected) % These variables must be initialised. Here or in the setupImpl function
        EKF;
        R;
        Q;
        gyro = [0;0;0];
        q;
        b_w = [0; 0;0];
        a_i;
        m_i;
    end
    
    methods (Access = protected)
        
        function resetImpl(obj)
            obj.gyro = [0;0;0]; 
        end
        
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.Q  = diag([obj.sigma_2_w',obj.sigma_2_b_w']); % contious-time process noise covariance matrix
            obj.R = diag([obj.sigma_2_a',obj.sigma_2_yaw']); % continous-time measurement noise covariance matrix
            x0 = [0;0;0;0;0;0];
            P0 = [1,0,0,0,0,0;
                   0,1,0,0,0,0;
                   0,0,1,0,0,0;
                   0,0,0,1,0,0;
                   0,0,0,0,1,0;
                   0,0,0,0,0,1];
            obj.EKF = EKF(obj.sampleTime, obj.Q, obj.R, P0, x0);
            obj.a_i = [0; 0; obj.g_mps2];
            obj.m_i = [1; 0; 0]; % ideal 
            obj.q = [1;0;0;0];
        end
        
        function attitude = stepImpl(obj,u, y)
            % Implement algorithm.
            obj.performEKF(u, y);
            q_hat = obj.attitude();
            q = Quaternion(q_hat);
            gamma_i_ahrsQuaternion = rad2deg(q.quatToEuler());
            attitude = gamma_i_ahrsQuaternion;
        end
    end
    
    methods
        
        % Constructor must be empty for matlab.System. In Matlab call
        % initFilter after the object was created. In simulink setupImpl()
        % will be called
        function obj = AHRSMEKF()
             % Support name-value pair arguments when constructing object
        end
        
        function initFilter(obj, sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma_2_w, sigma_2_a, sigma_2_b_w, sigma_2_yaw)
            % norm(q) must be one!
             x0 = [0;0;0;b_w0];
            obj.sigma_2_w = sigma_2_w;
            obj.sigma_2_a = sigma_2_a;
            obj.sigma_2_b_w = sigma_2_b_w;
            obj.sigma_2_yaw = sigma_2_yaw;
            obj.q = q0; % init attitude
            obj.b_w = b_w0; % init gyro bias
            obj.Q  = diag([sigma_2_w',sigma_2_b_w']); % contious-time process noise covariance matrix
            obj.R = diag([sigma_2_a',sigma_2_yaw']); % continous-time measurement noise covariance matrix
            obj.sampleTime = sampleTime;
            obj.EKF = EKF(sampleTime, obj.Q, obj.R, P0, x0);
            obj.delta_i = delta_i;
            obj.delta_d = delta_d;
            obj.g_mps2 = g_mps2;
            obj.a_i = [0; 0; g_mps2];
            obj.m_i = [1; 0; 0]; % ideal 
        end
        
        % estimate states
        function K = performEKF(obj, u, y)
            
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
            
            myInterpretation = 1;
            if myInterpretation
                [y_hat, H_lin] = obj.y_hat_and_H_lin(obj.EKF.x, u);
                % measurement step
                obj.EKF.correctorStep(y, y_hat, H_lin); % calculate the new states [hat(a)(+); b_what(+)]
            else
                for i=1:2
                    if i == 1
                        % this is without centripetal compensation
                        % accelerometer
                        y_part = y(1:3);
                        
                        q0_hat = obj.q(1);
                        q1_hat = obj.q(2);
                        q2_hat = obj.q(3);
                        q3_hat = obj.q(4);
                        g_mps2 = obj.g_mps2;
                        H_lin =  [0,	-g_mps2*(1-2*(q2_hat^2+q1_hat^2)),	2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat),	0,	0,	0;
                                g_mps2*(1-2*(q2_hat^2+q1_hat^2)),	0,	-2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat),	0,	0,	0;
                                -2*g_mps2*(q2_hat*q3_hat+q0_hat*q1_hat),	2*g_mps2*(q1_hat*q3_hat-q0_hat*q2_hat),	0,	0,	0,	0];
                            
                        y_hat = obj.Re_to_b(obj.q, [0;0;g_mps2]);
                        obj.EKF.R_k = obj.R(1:3,1:3)/ obj.sampleTime;
                        
                        y_hat = y_hat + H_lin*obj.EKF.x;
                    else
                        % magnetometer
                        y_part = y(4:6); 
                        
                        delta_i = obj.delta_i;
                        delta_d = obj.delta_d;
                        
                        H_lin = [0,	2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)-2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)+sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)),	-sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))-2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat),	0,	0,	0;
                                -2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q3_hat+q0_hat*q2_hat)-sin((pi*delta_i)/180)*(1-2*(q2_hat^2+q1_hat^2)),	0,	-cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))+2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)+2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat),	0,	0,	0;
                                sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q1_hat^2))+2*sin((pi*delta_i)/180)*(q2_hat*q3_hat+q0_hat*q1_hat)-2*cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q1_hat*q2_hat-q0_hat*q3_hat),	cos((pi*delta_d)/180)*cos((pi*delta_i)/180)*(1-2*(q3_hat^2+q2_hat^2))-2*sin((pi*delta_i)/180)*(q1_hat*q3_hat-q0_hat*q2_hat)-2*sin((pi*delta_d)/180)*cos((pi*delta_i)/180)*(q0_hat*q3_hat+q1_hat*q2_hat),	0,	0,	0,	0];
                        B = [cos((pi*delta_d)/180)*cos((pi*delta_i)/180);
                            -sin((pi*delta_d)/180)*cos((pi*delta_i)/180);
                            -sin((pi*delta_i)/180)];
                        y_hat = obj.Re_to_b(obj.q, B);
                        obj.EKF.R_k = obj.R(4:6, 4:6) / obj.sampleTime;
                        
                        y_hat = y_hat + H_lin*obj.EKF.x;
                    end
                    obj.EKF.correctorStep(y_part, y_hat, H_lin);
                end
            end
            
            obj.updateAttitude();
            % Wieso nicht hier den reset machen, sondern weiter oben?
            %obj.EKF.x = obj.EKF.x .* 0; % reset; this is wrong. works only, when all states are errors
            %obj.EKF.x = obj.EKF.x.*[0;0;0;1;1;1]; %reset only alpha
            
            K = obj.EKF.K;
        end
        
        % this function must be reimplemented in the derived class, because
        % the states are at different positions
        function updateAttitude(obj)
            % calculate corrected attitude
            % Komisch ist, dass wenn die Lage berechnet wird und dann
            % gleich resettet wird, funktionierts nicht (diesen Teil oben
            % vor dem reset machen, funktioniert nicht)
            alpha = obj.EKF.x(1:3);
            variant1 = 0;
            if variant1
                % variant 1 [2], eq. 102:
                alpha = alpha/2; % use this convention to calculate from alpha to q(alpha)
                qAlpha = [1; alpha];
                q = Quaternion(obj.q);
                obj.q = q.mult(qAlpha);
            else
                % variant 2 [3], eq. 22: (should be equivalent to variant 1)
                qAlpha = [2; alpha];
                % 'Dot notation on function call return value is not
                % allowed.'. (Simulink)
                q = Quaternion(obj.q);
                obj.q = q.mult(qAlpha);
                % in [3], eq.22 the multiplication is defined in the other
                % order, but there is the scalar also q3 and not q0!
            end
            % warning('AHRSMEKF: muss hier nochmals normalisiert 
            % werden? Ja bei den obigen Varianten schon (da delta a nicht 
            % ein unit quaternion ist). Wenn man aus [3] die Gleichung 21d 
            % verwenden würde, müsste man nicht mehr normalisieren');% 
            % normalize quaternion
            obj.q = obj.q / norm(obj.q); 
            obj.b_w = obj.EKF.x(4:6);
        end
        
        % return states
        function [a, b_w] = states(obj)
            a = obj.EKF.x(1:3);
            b_w = obj.EKF.x(4:6);
        end
        
        function q = attitude(obj)
            q = obj.q;
        end
        
        function b_w = gyroBias(obj)
            b_w = obj.b_w;
        end
        
        % calculate q_dot
        function q_dot = q_dot_f(obj, u)
          
            w_phi_m = u(1); % measured angular rates
            w_theta_m = u(2);
            w_psi_m = u(3);
            
            q_0 = obj.q(1);
            q_1 = obj.q(2);
            q_2 = obj.q(3);
            q_3 = obj.q(4);
            
            b_w_phi = obj.b_w(1);
            b_w_theta = obj.b_w(2);
            b_w_psi = obj.b_w(3);
            
            % 1/2 w x q
            q_dot = [(-q_2*(w_theta_m-b_w_theta)-q_3*(w_psi_m-b_w_psi)-q_1*(w_phi_m-b_w_phi))/2;
                    (-q_3*(w_theta_m-b_w_theta)+q_2*(w_psi_m-b_w_psi)+q_0*(w_phi_m-b_w_phi))/2;
                    (q_0*(w_theta_m-b_w_theta)-q_1*(w_psi_m-b_w_psi)+q_3*(w_phi_m-b_w_phi))/2;
                    (q_1*(w_theta_m-b_w_theta)+q_0*(w_psi_m-b_w_psi)-q_2*(w_phi_m-b_w_phi))/2];
        end
        
        % rotating vector v from earth to body with angle q
        function v = Re_to_b(obj, q, v)
            q_0 = q(1);
            q_1 = q(2);
            q_2 = q(3);
            q_3 = q(4);
            v = [1-2*(q_3^2+q_2^2),	2*(q_0*q_3+q_1*q_2),	2*(q_1*q_3-q_0*q_2);
                2*(q_1*q_2-q_0*q_3),	1-2*(q_3^2+q_1^2),	2*(q_2*q_3+q_0*q_1);
                2*(q_1*q_3+q_0*q_2),	2*(q_2*q_3-q_0*q_1),	1-2*(q_2^2+q_1^2)] * v;
        end

        % x_dot and G were calculated in maxima:
        % "KalmanfilterLagebestimmungEKF_BARTONSystemmatrizen/Quaternionen.wxm"
        % Chapter: Standard Multiplicative Extended Kalman Filter
        
        % calculate the derivative of the states
        function [x_dot, G_lin] = x_dot_and_G_lin(obj, x,w)
            
            w_phi_m = w(1);
            w_theta_m = w(2);
            w_psi_m = w(3);

            b_w_phi = x(4);
            b_w_theta = x(5);
            b_w_psi = x(6);
            
            G_lin = [0,	w_psi_m-b_w_psi,	b_w_theta-w_theta_m,	-1,	0,	0;
                    b_w_psi-w_psi_m,	0,	w_phi_m-b_w_phi,	0,	-1,	0;
                    w_theta_m-b_w_theta,	b_w_phi-w_phi_m,	0,	0,	0,	-1;
                    0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0];
                
            x_dot = G_lin * x;

        end
        
        function [y_hat, H_lin] = y_hat_and_H_lin(obj, x, u)
            
            delta_i = obj.delta_i;
            delta_d = obj.delta_d;
            g_mps2 = obj.g_mps2;

            % x: state vector 7x1
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