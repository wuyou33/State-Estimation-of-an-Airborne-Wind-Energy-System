classdef Quaternion < handle
    properties
        q0;
        q1;
        q2;
        q3;
        
        roll_est_deg;
        pitch_est_deg;
        yaw_est_deg;
    end
    
    methods
        function obj = Quaternion(q)
            obj.q0 = q(1);
            obj.q1 = q(2);
            obj.q2 = q(3);
            obj.q3 = q(4);
        end
        
        % this * p
        % implemented like in Maley2013
        function quat = mult(obj, q)
            
            p = [obj.q0;
                obj.q1;
                obj.q2;
                obj.q3];
            q =    [q(1);
                    q(2);
                    q(3);
                    q(4)];
            
            quat = [ p(1)*q(1) - p(2:4)'*q(2:4);
                    p(1)*q(2:4) + q(1)*p(2:4) + cross(p(2:4), q(2:4))];
        end
        
        function gamma = quatToEuler(obj) % [rad]
            q0 = obj.q0;
            q1 = obj.q1;
            q2 = obj.q2;
            q3 = obj.q3;
            % Problem: Quaternion [0.707,0,0.707,0] is in euler:
            % [0,90°,0] or [-180°, 90°, -180°] not unique
            % solved by the condition below
            C_ned2b  = [ 1-2*(q2^2+q3^2)    2*(q1*q2+q3*q0)     2*(q1*q3-q2*q0); ...
                         2*(q1*q2-q3*q0)    1-2*(q1^2+q3^2)     2*(q2*q3+q1*q0); ...
                         2*(q1*q3+q2*q0)    2*(q2*q3-q1*q0)     1-2*(q1^2+q2^2)];

            roll_est  = atan2(C_ned2b(2,3),C_ned2b(3,3));
            pitch_est = asin( -C_ned2b(1,3) );
            yaw_est   = atan2(C_ned2b(1,2),C_ned2b(1,1)); % -180 <= yaw <= 180
            if abs(C_ned2b(1,3)) > 1 - 1e-5
                % Pitch=+/-90deg case.  Underdetermined, so assume roll is zero,
                % and solve for pitch and yaw as follows:
                roll_est   = 0;
                pitch_est  = atan2( -C_ned2b(1,3), C_ned2b(3,3) );
                yaw_est    = atan2( -C_ned2b(2,1), C_ned2b(2,2) );
            end

            gamma = [roll_est;
                     pitch_est;
                     yaw_est ];
        end
        
        function gamma = debugAngleDeg(obj)
            roll = obj.roll_est_deg;
            pitch = obj.pitch_est_deg;
            yaw = obj.yaw_est_deg;
            
            gamma = [roll;
                     pitch;
                     yaw];
        end
    end
end