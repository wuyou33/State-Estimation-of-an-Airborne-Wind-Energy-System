classdef INSMEKF < AHRSMEKF
    %% Multiplicative EKF AHRS
    % Source: 
    % - [2]: Maley2013  - Multiplicative Quaternion Extended Kalman Filtering for Nonspinning Guided Projectiles
    % - [3]: Markley2003 - Attitude Error Representations for Kalman Filtering
    % - [4]: Patrick Spieler - https://github.com/Murmele/adcs
    % - [5]: Mattia Giurato - https://github.com/Murmele/Attitude_Estimation
    % 
    properties
    end
    
    methods
        
        function obj = INSMEKF()
            obj@AHRSMEKF();
        end
        
        function initFilter(obj, sampleTime, delta_i, delta_d, g_mps2, q0, b_w0, x0, P0, sigma2_process_noise, sigma2_measurement_noise)
           
            obj.q = q0; % init attitude
            obj.b_w = b_w0; % init gyro bias
            obj.Q  = sigma2_process_noise;
            obj.R = sigma2_measurement_noise;
            obj.sampleTime = sampleTime;
            obj.EKF = EKF(sampleTime, obj.Q, obj.R, P0, x0);
            obj.delta_i = delta_i;
            obj.delta_d = delta_d;
            obj.g_mps2 = g_mps2;
            obj.a_i = [0; 0; g_mps2];
            obj.m_i = [1; 0; 0]; % ideal 
        end
        
        function [gamma, pos, vel, b_w] = insStates(obj)
            q = Quaternion(obj.q);
            
            gamma = q.quatToEuler();
            pos = obj.EKF.x(4:6);
            vel = obj.EKF.x(7:9);
            b_w = obj.EKF.x(10:12);
            
        end
        
        % estimate states
        function K = performEKF(obj, u, y)
            
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
            
            [y_hat, H_lin] = obj.y_hat_and_H_lin(obj.EKF.x, u);
            % measurement step
            obj.EKF.correctorStep(y, y_hat, H_lin); % calculate the new states [hat(a)(+); b_what(+)]
            
            obj.updateAttitude();
            % Wieso nicht hier den reset machen, sondern weiter oben?
            %obj.EKF.x = obj.EKF.x .* 0; % reset; this is wrong. works only, when all states are errors
            %obj.EKF.x = obj.EKF.x.*[0;0;0;1;1;1;1;1;1;1;1;1]; %reset only alpha
            
            K = obj.EKF.K;
        end

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
            obj.b_w = obj.EKF.x(10:12);
        end

        
        % x_dot and G were calculated in maxima:
        % "EKFQuaternionen.wxm"
        % Chapter: Standard Multiplicative Extended Kalman Filter
        
        % calculate the derivative of the states
        function [x_dot, G_lin] = x_dot_and_G_lin(obj, x, u)
            
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
            a_bx = u(4);
            a_by = u(5);
            a_bz = u(6);
            
            q0_hat = obj.q(1);
            q1_hat = obj.q(2);
            q2_hat = obj.q(3);
            q3_hat = obj.q(4);
            
            g_mps2 = obj.g_mps2;
            
            x_dot = [a2*(b_w_theta-w_theta_m)+a1*(w_psi_m-b_w_psi)-b_w_phi;
                    a0*(b_w_psi-w_psi_m)+a2*(w_phi_m-b_w_phi)-b_w_theta;
                    a0*(w_theta_m-b_w_theta)+a1*(b_w_phi-w_phi_m)-b_w_psi;
                    vx;
                    vy;
                    vz;
                    (a1*a_bz-a2*a_by+a_bx)*(1-2*(q3_hat^2+q2_hat^2))+2*(a_bz+a0*a_by-a1*a_bx)*(q1_hat*q3_hat+q0_hat*q2_hat)+2*(-a0*a_bz+a_by+a2*a_bx)*(q1_hat*q2_hat-q0_hat*q3_hat);
                    (-a0*a_bz+a_by+a2*a_bx)*(1-2*(q3_hat^2+q1_hat^2))+2*(a_bz+a0*a_by-a1*a_bx)*(q2_hat*q3_hat-q0_hat*q1_hat)+2*(a1*a_bz-a2*a_by+a_bx)*(q0_hat*q3_hat+q1_hat*q2_hat);
                    2*(-a0*a_bz+a_by+a2*a_bx)*(q2_hat*q3_hat+q0_hat*q1_hat)+2*(a1*a_bz-a2*a_by+a_bx)*(q1_hat*q3_hat-q0_hat*q2_hat)+(a_bz+a0*a_by-a1*a_bx)*(1-2*(q2_hat^2+q1_hat^2))-g_mps2;
                    0;
                    0;
                    0];
            G_lin = [0,	w_psi_m-b_w_psi,	b_w_theta-w_theta_m,	0,	0,	0,	0,	0,	0,	-1,	a2,	-a1;
                    b_w_psi-w_psi_m,	0,	w_phi_m-b_w_phi,	0,	0,	0,	0,	0,	0,	-a2,	-1,	a0;
                    w_theta_m-b_w_theta,	b_w_phi-w_phi_m,	0,	0,	0,	0,	0,	0,	0,	a1,	-a0,	-1;
                    0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0;
                    2*a_by*(q1_hat*q3_hat+q0_hat*q2_hat)-2*a_bz*(q1_hat*q2_hat-q0_hat*q3_hat),	a_bz*(1-2*(q3_hat^2+q2_hat^2))-2*a_bx*(q1_hat*q3_hat+q0_hat*q2_hat),	2*a_bx*(q1_hat*q2_hat-q0_hat*q3_hat)-a_by*(1-2*(q3_hat^2+q2_hat^2)),	0,	0,	0,	0,	0,	0,	0,	0,	0;
                    2*a_by*(q2_hat*q3_hat-q0_hat*q1_hat)-a_bz*(1-2*(q3_hat^2+q1_hat^2)),	2*a_bz*(q0_hat*q3_hat+q1_hat*q2_hat)-2*a_bx*(q2_hat*q3_hat-q0_hat*q1_hat),	a_bx*(1-2*(q3_hat^2+q1_hat^2))-2*a_by*(q0_hat*q3_hat+q1_hat*q2_hat),	0,	0,	0,	0,	0,	0,	0,	0,	0;
                    a_by*(1-2*(q2_hat^2+q1_hat^2))-2*a_bz*(q2_hat*q3_hat+q0_hat*q1_hat),	2*a_bz*(q1_hat*q3_hat-q0_hat*q2_hat)-a_bx*(1-2*(q2_hat^2+q1_hat^2)),	2*a_bx*(q2_hat*q3_hat+q0_hat*q1_hat)-2*a_by*(q1_hat*q3_hat-q0_hat*q2_hat),	0,	0,	0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0;
                    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0];
        end
        
        function [y_hat, H_lin] = y_hat_and_H_lin(obj, x, u)
            
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