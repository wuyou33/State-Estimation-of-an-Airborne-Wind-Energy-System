% Returns a coord system which can be plotted. Rotate and translate coord
% system.
% Definition: North-East-Down (NED)

% test:
% S = coordSystem(1);
% S.edit([1;0;0], [90*pi/180;0;0]);

classdef CoordSystem < handle
    properties
        achsX;
        achsY;
        achsZ;
        length;
    end
    
    methods
        function this = CoordSystem(arrowLength)
            this.length = arrowLength;
        end
        
        function [achsX, achsY, achsZ] = editDeg(obj, pos, gamma)
            gamma = gamma * pi / 180;
            [achsX, achsY, achsZ] = edit(obj, pos, gamma);
        end
        
        function [achsX, achsY, achsZ] = edit(obj, pos, gamma)
            e1 = obj.length * [1;0;0]; % north
            e2 = obj.length * [0;1;0]; % east
            e3 = obj.length * [0;0;-1]; % down
            pos_e(:,1) = obj.transformationMatrix(pos, gamma)*[e1;1];
            pos_e(:,2) = obj.transformationMatrix(pos, gamma)*[e2;1];
            pos_e(:,3) = obj.transformationMatrix(pos, gamma)*[e3;1];
            pos_e = pos_e(1:3,:);
            % equivalent
            %pos_e2 = obj.rotationMatrix(gamma)*[e1, e2, e3]+ [pos, pos, pos];
            
            pos_e_original = [pos + e1, pos + e2, pos + e3];
            
            pos_origX = [pos, pos_e_original(:, 1)];
            pos_origY = [pos, pos_e_original(:, 2)];
            pos_origZ = [pos, pos_e_original(:, 3)];
            
            obj.achsX = [pos, pos_e(:, 1)];
            obj.achsY = [pos, pos_e(:, 2)];
            obj.achsZ = [pos, pos_e(:, 3)];
            
            achsX = obj.achsX;
            achsY = obj.achsY;
            achsZ = obj.achsZ;
            
%             % for testing
%             f = figure();
%             xlabel('X');
%             ylabel('Y');
%             zlabel('Z');
%             hold('on');
%             plot3(pos_origX(1,:), pos_origX(2,:), pos_origX(3,:),'m', 'DisplayName', 'Xorig');
%             plot3(pos_origY(1,:), pos_origY(2,:), pos_origY(3,:),'y', 'DisplayName', 'Yorig');
%             plot3(pos_origZ(1,:), pos_origZ(2,:), pos_origZ(3,:),'c', 'DisplayName', 'Zorig');
%             plot3(obj.achsX(1,:), obj.achsX(2,:), obj.achsX(3,:),'r', 'DisplayName', 'X');
%             plot3(obj.achsY(1,:), obj.achsY(2,:), obj.achsY(3,:),'g', 'DisplayName', 'Y');
%             plot3(obj.achsZ(1,:), obj.achsZ(2,:), obj.achsZ(3,:),'b', 'DisplayName', 'Z');
%             hold('off');
        end
        
        function [achsX, achsY, achsZ] = axis(obj);
            achsX = obj.achsX;
            achsY = obj.achsY;
            achsZ = obj.achsZ;
        end
        
        function T = transformationMatrix(obj, pos, gamma)
            T = [obj.rotationMatrix(gamma), pos;
                zeros(1,3), 1];
        end
        
        function M = rotationMatrix(obj, gamma)
            phi = gamma(1);
            theta = gamma(2);
            psi = gamma(3);

            % Rb_to_e
            M = [cos(psi)*cos(theta),	sin(phi)*cos(psi)*sin(theta)-cos(phi)*sin(psi),	cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi);
                sin(psi)*cos(theta),	sin(phi)*sin(psi)*sin(theta)+cos(phi)*cos(psi),	cos(phi)*sin(psi)*sin(theta)-sin(phi)*cos(psi);
                -sin(theta),	sin(phi)*cos(theta),	cos(phi)*cos(theta)];
        end
    end
end