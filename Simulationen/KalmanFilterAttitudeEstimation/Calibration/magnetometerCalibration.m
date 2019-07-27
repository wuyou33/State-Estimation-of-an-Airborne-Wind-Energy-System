% source:
% - [1]: 838300 - Applications of magnetic sensors for low cost compass systems
% - [2]: Ozyagcilar2015a - Implementing a Tilt-Compensated eCompass using Accelerometer and Magnetometer Sensors
% - [3]: Ozyagcilar2015 - Calibrating an eCompass in the Presence of Hard- and Soft-Iron Interference

% - mag_m_unnormalized: 3x1 unnormalized magnetometer values
% - plotting: if calibration should be plotted or not (for demonstration
% purpose only)
% hardIronCalibrationMethod: different methods to determine the hardIron
% Calibration {freescale (fast and accurate), ellipsoid_fit, simple (default, simple, less accurate)}
% softIronCalibrationMethod: different methods to determine the
% softIronCalibration {freescale (not working), ownInterpretation (my
% interpretation from different sources, default)}

% offset: hard iron calibration: mag_m = mag_m_meas - offset (mag_m are the
% calibrated values, mag_m_meas are the uncalibrated values)
% softIronMatrix: mag_m = softIronMatrix * (mag_m_meas - offset)
function [offset, softIronMatrix] = magnetometerCalibration(mag_m_unnormalized, plotting, hardIronCalibrationMethod, softIronCalibrationMethod)

    circle_size = 5; % used in the scatter plot
    circle_color = 'filled'; % used in the scatter plot
    
    % circles in the plot ( not used at the moment)
    folderNameFlatRotation = {}; %{'YawRotation_Pitch0_Roll0', 'YawRotation_Pitch0_Roll90', 'YawRotation_Pitch-90_Roll0'}; % name of the folder which data should be used
    labelsFlatRotation = {'Pitch=0°, Roll=0°', 'Pitch=0°, Roll=90°', 'Pitch=-90°, Roll=0°'};
    style = {'r', 'g', 'b'};

	if strcmp(softIronCalibrationMethod, '')
		softIronCalibrationMethod = 'ownInterpretation';
    end
    
    if strcmp(softIronCalibrationMethod, 'freescale') == 0 && strcmp(softIronCalibrationMethod, 'ownInterpretation') == 0
        error(['no correction calibration Method selected: ', softIronCalibrationMethod]);
    end

    %% Hard iron compensation ( All three methods should be equivalent, but second method is not so accurate, but the fastest)
    
        %% first method: Freescale/NXP Paper [3]
        X = [mag_m_unnormalized', ones(length(mag_m_unnormalized), 1)]; % [3], eq. 30
        Y = (vecnorm(mag_m_unnormalized).^2)'; % [3], eq. 29
        beta = (X'*X)\(X'*Y); % [3], eq. 31
        % better than using the mean, because this is least square estimation
        % same result than ellipsoid fit below
        offset_freescale = beta(1:3)/2;
        % magnetic field strength
        B_square = beta(4) + offset_freescale(1)^2 + offset_freescale(2)^2 + offset_freescale(3)^2;
        
        %% second method: calculating offset from min and max values
        maxX = max(mag_m_unnormalized(1,:));
        minX = min(mag_m_unnormalized(1,:));
        maxY = max(mag_m_unnormalized(2,:));
        minY = min(mag_m_unnormalized(2,:));
        maxZ = max(mag_m_unnormalized(3,:));
        minZ = min(mag_m_unnormalized(3,:));

        % easiest method to calculate the hard iron offset.
        % Other methods are the method from freescale [3]. See implementation
        % above or directly from the elipsoid fitting below
        offset_simple(1,:) = (maxX + minX) / 2; % Offset X
        offset_simple(2,:) = (maxY + minY) / 2; % Offset Y
        offset_simple(3,:) = (maxZ + minZ) / 2; % Offset Z
    
        %% third method: calculate offset from ellipsoid_fit
        % ellipsoid_fit (needed later also in the soft iron compensation
        [ center, radii, evecs, evals, v, chi2 ] = ellipsoid_fit( mag_m_unnormalized', '' );
    
    
    % different methods to determine the hard iron offset
    %hardIronCalibrationMethod = 'freescale'; % 'freescale', 'ellipsoid_fit', 'simple'
    display(['Use "' hardIronCalibrationMethod '" method to determine hard iron offset']);
    if (strcmp(hardIronCalibrationMethod, 'freescale'))
        offset = offset_freescale;
    elseif (strcmp(hardIronCalibrationMethod, 'ellipsoid_fit'))
        offset = center;
    else
        offset = offset_simple; % default if nothing selected
    end
    
    clear maxX minX maxY minY maxZ minZ
    

    %% Soft iron compensation (exist also multiple solutions to calculate the soft iron matrix)
    
    
    %% First method: from [3]
    % Eq. 11 solved to A = (W^-1)^T W^-1
    % finding A from the above ellipsoid_fit with the eigenvalues and the
    % eigenvectors
    % https://math.stackexchange.com/questions/54818/construct-matrix-given-eigenvalues-and-eigenvectors
    %A*V = V*D % see matlab eig
    if strcmp(softIronCalibrationMethod, 'freescale')
        D = evals; % eigenvalue matrix
        V = evecs;
        A = V*D*inv(V);
        softIronMatrix = abs(sqrt(A)); % hat A eventuell negative einträge? Dürfte nicht sein, da sonst die Wurzel daraus zu Imaginären Werten führen würde
        mag_calib = softIronMatrix * ((mag_m_unnormalized - offset));
        warning(['Seems not to be correct,' ...
                'the results are ellipsoids and no sphere!']);
    end
    % softIronMatrix = W_1;
    
    %% Second Method: own interpretation:
    % eigenvectors are the main axis of the ellipsoid
    % 1) rotate ellipsoid to match the axis x,y,z
    % 2) scale ellipsoid with the radius
    % 3) rotate ellipsoid back to match the previous axis
    % 4) rotate about the circle center lines, so the circles are mapped to
    % the planes (front, up, right)
    % (this is only correct if roll and pitch angles are zero!)

    
    mag_m_hard_iron_comp = mag_m_unnormalized - offset; % important to subtract offset, otherwise the eigenvectors are plotted wrong
    if strcmp(softIronCalibrationMethod, 'ownInterpretation')
        
        f_ownInterpretation = figure('Name', 'Magnetometer calibration, own interpretation');
        figure(f_ownInterpretation);
        movegui(f_ownInterpretation, 'south');
        subplot(2,3,1);
        pbaspect([1 1 1]);
        xlabel('Xh');
        ylabel('Yh');
        zlabel('Zh');
        title('Without rotation');
        hold on;
        scatter3(mag_m_hard_iron_comp(1,:), mag_m_hard_iron_comp(2, :), mag_m_hard_iron_comp(3, :), circle_size, circle_color, 'DisplayName', 'Raw data');
        % draw eigenvectors:
        [ center, radii, evecs, evals, v, chi2 ] = ellipsoid_fit( mag_m_hard_iron_comp', '' );
        axis1 = [center, evecs(:, 1)*radii(1)];
        axis2 = [center, evecs(:, 2)*radii(2)];
        axis3 = [center, evecs(:, 3)*radii(3)];
        plot3(axis1(1, :), axis1(2, :), axis1(3, :), 'r','HandleVisibility','off');
        plot3(axis2(1, :), axis2(2, :), axis2(3, :), 'g','HandleVisibility','off');
        plot3(axis3(1, :), axis3(2, :), axis3(3, :), 'b','HandleVisibility','off');
        for i=1:length(folderNameFlatRotation)
            sensorData = readArduinoData(fullfile(arduinoPath, folderNameFlatRotation{i}), 'config.m', 'received.log', g_mps2);
            magData = sensorData.meas.m_b_unnormalized - offset;
            plot3(magData(1, :), magData(2, :), magData(3, :), style{i}, 'DisplayName', labelsFlatRotation{i});
            [circleCenter, CircleRadius, CircleNormalvector] =  findCircleIn3DLeastSquare(magData);
            centerVector = [offset, circleCenter];
            plot3(centerVector(1, :),centerVector(2, :),centerVector(3, :) , style{i},'HandleVisibility','off');
        end
        hold off;
        ax = gca; % current axes
        ax_lim_max = max([abs(ax.XLim), abs(ax.YLim), abs(ax.ZLim)]);
        ax.XLim = [-ax_lim_max, ax_lim_max];
        ax.YLim = [-ax_lim_max, ax_lim_max];
        ax.ZLim = [-ax_lim_max, ax_lim_max];
        
        % 1)
        subplot(2,3,2);
        pbaspect([1,1,1]);
        title('With rotation');
        xlabel('Xh');
        ylabel('Yh');
        zlabel('Zh');
        hold on;
        R_ev = evecs;
        mag_m_unnormalized_rot = R_ev'*mag_m_hard_iron_comp;
        scatter3(mag_m_unnormalized_rot(1,:), mag_m_unnormalized_rot(2, :), mag_m_unnormalized_rot(3, :), circle_size, circle_color, 'DisplayName', 'Raw data');
        % draw eigenvectors:
        %[ center, radii, evecs, evals, v, chi2 ] = ellipsoid_fit( mag_m_unnormalized_rot', '' );
        evecs_new = R_ev'*evecs;
        axis1 = [center, evecs_new(:, 1)*radii(1)];
        axis2 = [center, evecs_new(:, 2)*radii(2)];
        axis3 = [center, evecs_new(:, 3)*radii(3)];
        plot3(axis1(1, :), axis1(2, :), axis1(3, :), 'r','HandleVisibility','off');
        plot3(axis2(1, :), axis2(2, :), axis2(3, :), 'g','HandleVisibility','off');
        plot3(axis3(1, :), axis3(2, :), axis3(3, :), 'b','HandleVisibility','off');
        for i=1:length(folderNameFlatRotation)
            sensorData = readArduinoData(fullfile(arduinoPath, folderNameFlatRotation{i}), 'config.m', 'received.log', g_mps2);
            magData = R_ev'*(sensorData.meas.m_b_unnormalized - offset);
            plot3(magData(1, :), magData(2, :), magData(3, :), style{i}, 'DisplayName', labelsFlatRotation{i});
            [circleCenter, CircleRadius, CircleNormalvector] =  findCircleIn3DLeastSquare(magData);
            centerVector = [offset, circleCenter];
            plot3(centerVector(1, :),centerVector(2, :),centerVector(3, :) , style{i},'HandleVisibility','off');
        end
        hold off;
        ax = gca; % current axes
        ax_lim_max = max([abs(ax.XLim), abs(ax.YLim), abs(ax.ZLim)]);
        ax.XLim = [-ax_lim_max, ax_lim_max];
        ax.YLim = [-ax_lim_max, ax_lim_max];
        ax.ZLim = [-ax_lim_max, ax_lim_max];
        
        % 2)
        subplot(2,3,3);
        pbaspect([1,1,1]);
        title('With rotation and scale');
        xlabel('Xh');
        ylabel('Yh');
        zlabel('Zh');
        hold on;
        R_ev = evecs;
        radiiMatrixInv = diag(1./radii);
        mag_m_unnormalized_rot_scale = radiiMatrixInv*mag_m_unnormalized_rot;
        scatter3(mag_m_unnormalized_rot_scale(1,:), mag_m_unnormalized_rot_scale(2, :), mag_m_unnormalized_rot_scale(3, :), circle_size, circle_color, 'DisplayName', 'Raw data');
        % draw eigenvectors:
        %[ center, radii, evecs, evals, v, chi2 ] = ellipsoid_fit( mag_m_unnormalized_rot', '' );
        evecs_new = R_ev'*evecs;
        centerZero = [0;0;0];
        axis1 = [centerZero, evecs_new(:, 1)];
        axis2 = [centerZero, evecs_new(:, 2)];
        axis3 = [centerZero, evecs_new(:, 3)];
        plot3(axis1(1, :), axis1(2, :), axis1(3, :), 'r','HandleVisibility','off');
        plot3(axis2(1, :), axis2(2, :), axis2(3, :), 'g','HandleVisibility','off');
        plot3(axis3(1, :), axis3(2, :), axis3(3, :), 'b','HandleVisibility','off');
        for i=1:length(folderNameFlatRotation)
            sensorData = readArduinoData(fullfile(arduinoPath, folderNameFlatRotation{i}), 'config.m', 'received.log', g_mps2);
            magData = radiiMatrixInv*R_ev'*(sensorData.meas.m_b_unnormalized - offset);
            plot3(magData(1, :), magData(2, :), magData(3, :), style{i}, 'DisplayName', labelsFlatRotation{i});
            [circleCenter, CircleRadius, CircleNormalvector] =  findCircleIn3DLeastSquare(magData);
            centerVector = [[0;0;0], circleCenter];
            plot3(centerVector(1, :),centerVector(2, :),centerVector(3, :) , style{i},'HandleVisibility','off');
        end
        hold off;
        ax = gca; % current axes
        ax_lim_max = max([abs(ax.XLim), abs(ax.YLim), abs(ax.ZLim)]);
        ax.XLim = [-ax_lim_max, ax_lim_max];
        ax.YLim = [-ax_lim_max, ax_lim_max];
        ax.ZLim = [-ax_lim_max, ax_lim_max];
        
        % 3)
        subplot(2,3,4);
        pbaspect([1 1 1]);
        xlabel('Xh');
        ylabel('Yh');
        zlabel('Zh');
        title('Scaled matrix rotated back');
        hold on;
        % rotate back
        mag_m_unnormalized_rot_scale_rot = R_ev*mag_m_unnormalized_rot_scale;
        scatter3(mag_m_unnormalized_rot_scale_rot(1,:), mag_m_unnormalized_rot_scale_rot(2, :), mag_m_unnormalized_rot_scale_rot(3, :), circle_size, circle_color, 'DisplayName', 'Raw data');
        evecs_new_rot_back = R_ev*evecs_new;
        centerZero = [0;0;0];
        axis1 = [centerZero, evecs_new_rot_back(:, 1)];
        axis2 = [centerZero, evecs_new_rot_back(:, 2)];
        axis3 = [centerZero, evecs_new_rot_back(:, 3)];
        plot3(axis1(1, :), axis1(2, :), axis1(3, :), 'r','HandleVisibility','off');
        plot3(axis2(1, :), axis2(2, :), axis2(3, :), 'g','HandleVisibility','off');
        plot3(axis3(1, :), axis3(2, :), axis3(3, :), 'b','HandleVisibility','off');
        for i=1:length(folderNameFlatRotation)
            sensorData = readArduinoData(fullfile(arduinoPath, folderNameFlatRotation{i}), 'config.m', 'received.log', g_mps2);
            magData = R_ev*radiiMatrixInv*R_ev'*(sensorData.meas.m_b_unnormalized - offset);
            plot3(magData(1, :), magData(2, :), magData(3, :), style{i}, 'DisplayName', labelsFlatRotation{i});
            [circleCenter, CircleRadius, CircleNormalvector] =  findCircleIn3DLeastSquare(magData);
            centerVector = [[0;0;0], circleCenter];
            plot3(centerVector(1, :),centerVector(2, :),centerVector(3, :) , style{i},'HandleVisibility','off');
        end
        hold off;
        ax = gca; % current axes
        ax_lim_max = max([abs(ax.XLim), abs(ax.YLim), abs(ax.ZLim)]);
        ax.XLim = [-ax_lim_max, ax_lim_max];
        ax.YLim = [-ax_lim_max, ax_lim_max];
        ax.ZLim = [-ax_lim_max, ax_lim_max];
        
        % resulting matrix when combining all above steps
        softIronMatrix = R_ev*radiiMatrixInv*R_ev';
        
        % 4) (seems not to be correct)
        % maybe this is, because the circles must not be completely aligned
        % with the axis, because while rotaing the board was tilded a
        % little bit
%         arduinoPath = '../../Arduino10DOFSensor/DataAndMatlabImport/OutsideNewComplete';
%         addpath(arduinoPath);
%         subplot(2,3,5);
%         pbaspect([1 1 1]);
%         xlabel('Xh');
%         ylabel('Yh');
%         zlabel('Zh');
%         title('Scaled matrix rotated back and rotated, that circles match');
%         hold on;
%         
%         % get centers of the circles
%         rotationsCenter = zeros(3, 3);
%         for i=1:length(folderNameFlatRotation)
%             sensorData = readArduinoData(fullfile(arduinoPath, folderNameFlatRotation{i}), 'config.m', 'received.log', g_mps2);
%             magData = R_ev*radiiMatrixInv*R_ev'*(sensorData.meas.m_b_unnormalized - offset);
%             [circleCenter, CircleRadius, CircleNormalvector] =  findCircleIn3DLeastSquare(magData);
%             rotationsCenter(:, 3-i+1) = circleCenter;
%         end
%         mag_m_unnormalized_rot_scale_rot = rotationsCenter'*R_ev*mag_m_unnormalized_rot_scale;
%         scatter3(mag_m_unnormalized_rot_scale_rot(1,:), mag_m_unnormalized_rot_scale_rot(2, :), mag_m_unnormalized_rot_scale_rot(3, :), circle_size, circle_color, 'DisplayName', 'Raw data');
%         evecs_new_rot_back = rotationsCenter'*R_ev*evecs_new;
%         centerZero = [0;0;0];
%         axis1 = [centerZero, evecs_new_rot_back(:, 1)];
%         axis2 = [centerZero, evecs_new_rot_back(:, 2)];
%         axis3 = [centerZero, evecs_new_rot_back(:, 3)];
%         plot3(axis1(1, :), axis1(2, :), axis1(3, :), 'r','HandleVisibility','off');
%         plot3(axis2(1, :), axis2(2, :), axis2(3, :), 'g','HandleVisibility','off');
%         plot3(axis3(1, :), axis3(2, :), axis3(3, :), 'b','HandleVisibility','off');
%         for i=1:length(folderNameFlatRotation)
%             sensorData = readArduinoData(fullfile(arduinoPath, folderNameFlatRotation{i}), 'config.m', 'received.log', g_mps2);
%             magData = rotationsCenter'*R_ev*radiiMatrixInv*R_ev'*(sensorData.meas.m_b_unnormalized - offset);
%             plot3(magData(1, :), magData(2, :), magData(3, :), style{i}, 'DisplayName', labelsFlatRotation{i});
%             [circleCenter, CircleRadius, CircleNormalvector] =  findCircleIn3DLeastSquare(magData);
%             centerVector = [[0;0;0], circleCenter];
%             plot3(centerVector(1, :),centerVector(2, :),centerVector(3, :) , style{i},'HandleVisibility','off');
%         end
%         hold off;
%         ax = gca; % current axes
%         ax_lim_max = max([abs(ax.XLim), abs(ax.YLim), abs(ax.ZLim)]);
%         ax.XLim = [-ax_lim_max, ax_lim_max];
%         ax.YLim = [-ax_lim_max, ax_lim_max];
%         ax.ZLim = [-ax_lim_max, ax_lim_max];
    end
    %% Plotting 
    %warning('Plotting is off, for debugging. Remove this set to zero!');
    %plotting = 0;
    if plotting
        
        %% Plotting uncompensated / Hard iron compensated Magnetometer values
        f_uncomp = figure('Name', 'Magnetometer 2D');
        figure(f_uncomp);
        movegui(f_uncomp,'northwest');
        subplot(2,2,1);
        pbaspect([1 1 1])
        hold on;
        scatter3(mag_m_unnormalized(1,:), mag_m_unnormalized(2, :), mag_m_unnormalized(3, :), circle_size, circle_color, 'DisplayName', 'Raw data');
        
        % source: https://diydrones.com/forum/topics/magnetometer-soft-and-hard-iron-calibration?commentId=705844%3AComment%3A786369
        % The grey line is arbitrary rotations and the red, green and blue lines are rotations on a flat level surface around x, y and z (axis pointing down each time).

        display('Die Kreise sind nicht exakt in der ebene (Dies stellt den Softiron fehler dar! Aber auch nur, wenn perfekt in der Ebene gedreht wurde!)');
        arduinoPath = '../../Arduino10DOFSensor/DataAndMatlabImport/OutsideNewComplete';
        addpath(arduinoPath);

        for i=1:length(folderNameFlatRotation)
            sensorData = readArduinoData(fullfile(arduinoPath, folderNameFlatRotation{i}), 'config.m', 'received.log', g_mps2);
            magData = sensorData.meas.m_b_unnormalized;
            plot3(magData(1, :), magData(2, :), magData(3, :), style{i}, 'DisplayName', labelsFlatRotation{i});
            [circleCenter, CircleRadius, CircleNormalvector] =  findCircleIn3DLeastSquare(magData);
            centerVector = [offset, circleCenter];
            plot3(centerVector(1, :),centerVector(2, :),centerVector(3, :) , style{i},'HandleVisibility','off');
        end

        clear magData circleCenter radius normalvector centerVector
        % make that all axis have same limit
        ax = gca; % current axes
        ax_lim_max = max([abs(ax.XLim), abs(ax.YLim), abs(ax.ZLim)]);
        ax.XLim = [-ax_lim_max, ax_lim_max];
        ax.YLim = [-ax_lim_max, ax_lim_max];
        ax.ZLim = [-ax_lim_max, ax_lim_max];
        xlabel('Xh');
        ylabel('Yh');
        zlabel('Zh');
        hold off;
        subplot(2,2,2);
        pbaspect([1 1 1])
        hold on;
        scatter(mag_m_unnormalized(1,:), mag_m_unnormalized(2,:), circle_size, circle_color, 'DisplayName', 'Raw data');
        scatter(mag_m_unnormalized(1,:) - offset(1), mag_m_unnormalized(2,:) - offset(2), circle_size, circle_color,'DisplayName', 'With Hard iron compensation');
        % make that all axis have same limit
        ax = gca; % current axes
        ax_lim_max = max([abs(ax.XLim), abs(ax.YLim)]);
        ax.XLim = [-ax_lim_max, ax_lim_max];
        ax.YLim = [-ax_lim_max, ax_lim_max];
        
        hold off;
        title('XY-Plane (Top/Bottom)');
        xlabel('Xh');
        ylabel('Yh');
        subplot(2,2,3);
        pbaspect([1 1 1])
        hold on;
        scatter(mag_m_unnormalized(1,:), mag_m_unnormalized(3,:), circle_size, circle_color, 'DisplayName', 'Raw data');
        scatter(mag_m_unnormalized(1,:) - offset(1), mag_m_unnormalized(3,:) - offset(3), circle_size, circle_color, 'DisplayName', 'With Hard iron compensation');
        % make that all axis have same limit
        ax = gca; % current axes
        ax_lim_max = max([abs(ax.XLim), abs(ax.YLim)]);
        ax.XLim = [-ax_lim_max, ax_lim_max];
        ax.YLim = [-ax_lim_max, ax_lim_max];
        hold off;
        title('XZ-Plane (Front)');
        xlabel('Xh');
        ylabel('Zh');
        subplot(2,2,4);
        title('YZ-Plane (Side)');
        xlabel('Yh');
        ylabel('Zh');
        pbaspect([1 1 1])
        hold on;
        scatter(mag_m_unnormalized(2,:), mag_m_unnormalized(3,:), circle_size, circle_color, 'DisplayName', 'Raw data');
        scatter(mag_m_unnormalized(2,:) - offset(2), mag_m_unnormalized(3,:) - offset(3), circle_size, circle_color, 'DisplayName', 'With Hard iron compensation');
        % make that all axis have same limit
        ax = gca; % current axes
        ax_lim_max = max([abs(ax.XLim), abs(ax.YLim)]);
        ax.XLim = [-ax_lim_max, ax_lim_max];
        ax.YLim = [-ax_lim_max, ax_lim_max];
        hold off;
        
        % plot eigenvectors and isosurface of the ellipsoid
        subplot(2,2,1);
        hold on;
        mind = min( mag_m_unnormalized' );
        maxd = max( mag_m_unnormalized' );
        nsteps = 50;
        step = ( maxd - mind ) / nsteps;
        [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
        Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
              2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
              2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
        %p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );
        %set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );

        % draw eigenvectors:
        axis1 = [center, evecs(:, 1)*radii(1)];
        axis2 = [center, evecs(:, 2)*radii(2)];
        axis3 = [center, evecs(:, 3)*radii(3)];
        plot3(axis1(1, :), axis1(2, :), axis1(3, :), 'r','HandleVisibility','off');
        plot3(axis2(1, :), axis2(2, :), axis2(3, :), 'g','HandleVisibility','off');
        plot3(axis3(1, :), axis3(2, :), axis3(3, :), 'b','HandleVisibility','off');
        hold off;

        %% plot after soft and hard iron calibration
%         f_mag_comp = figure('Name', 'Magnetometer after soft/hard iron calibration');
%         figure(f_mag_comp);
%         movegui(f_mag_comp,'southwest');
%         subplot(2,2,1);
%         pbaspect([1 1 1])
%         hold on;
%         scatter3(mag_calib(1,:), mag_calib(2, :), mag_calib(3, :), circle_size, circle_color);
%         % source: https://diydrones.com/forum/topics/magnetometer-soft-and-hard-iron-calibration?commentId=705844%3AComment%3A786369
%         % The grey line is arbitrary rotations and the red, green and blue lines are rotations on a flat level surface around x, y and z (axis pointing down each time).
%         arduinoPath = '../../Arduino10DOFSensor/DataAndMatlabImport/OutsideNewComplete';
%         addpath(arduinoPath);
%         
%         for i=1: length(folderNameFlatRotation)
%             sensorData = readArduinoData(fullfile(arduinoPath, folderNameFlatRotation{i}), 'config.m', 'received.log', g_mps2);
%             magData = sensorData.meas.m_b_unnormalized;
%             mag_comp = W_1*(magData-offset);
%             plot3(mag_comp(1, :), mag_comp(2, :), mag_comp(3, :), style{i},'DisplayName', labelsFlatRotation{i});
%         end
% 
%         clear mag_comp magData
%         
%         % draw eigenvectors:
%         axis1 = [zeros(3,1), evecs(:, 1)/norm(evecs(:,1))];
%         axis2 = [zeros(3,1), evecs(:, 2)/norm(evecs(:,2))];
%         axis3 = [zeros(3,1), evecs(:, 3)/norm(evecs(:,3))];
%         plot3(axis1(1, :), axis1(2, :), axis1(3, :), 'r');
%         plot3(axis2(1, :), axis2(2, :), axis2(3, :), 'g');
%         plot3(axis3(1, :), axis3(2, :), axis3(3, :), 'b');
%         % make that all axis have same limit
%         ax = gca; % current axes
%         ax_lim_max = max([abs(ax.XLim), abs(ax.YLim), abs(ax.ZLim)]);
%         ax.XLim = [-ax_lim_max, ax_lim_max];
%         ax.YLim = [-ax_lim_max, ax_lim_max];
%         ax.ZLim = [-ax_lim_max, ax_lim_max];
%         xlabel('Xh');
%         ylabel('Yh');
%         zlabel('Zh');
%         hold off;
%         subplot(2,2,2);
%         pbaspect([1 1 1])
%         title('XY-Plane (Top/Bottom)');
%         xlabel('Xh');
%         ylabel('Yh');
%         hold on;
%         scatter(mag_calib(1,:), mag_calib(2,:), circle_size, circle_color, 'DisplayName', 'soft/hard iron compensation');
%         plot(axis1(1, :), axis1(2, :), 'r','HandleVisibility','off');
%         plot(axis2(1, :), axis2(2, :), 'g','HandleVisibility','off');
%         plot(axis3(1, :), axis3(2, :), 'b','HandleVisibility','off');
%         % make that all axis have same limit
%         ax = gca; % current axes
%         ax_lim_max = max([abs(ax.XLim), abs(ax.YLim)]);
%         ax.XLim = [-ax_lim_max, ax_lim_max];
%         ax.YLim = [-ax_lim_max, ax_lim_max];
%         hold off;
%         subplot(2,2,3);
%         pbaspect([1 1 1])
%         title('XZ-Plane (Front)');
%         xlabel('Xh');
%         ylabel('Zh');
%         hold on;
%         scatter(mag_calib(1,:), mag_calib(3,:), circle_size, circle_color, 'DisplayName', 'soft/hard iron compensation');
%         plot(axis1(1, :), axis1(3, :), 'r','HandleVisibility','off');
%         plot(axis2(1, :), axis2(3, :), 'g','HandleVisibility','off');
%         plot(axis3(1, :), axis3(3, :), 'b','HandleVisibility','off');
%         % make that all axis have same limit
%         ax = gca; % current axes
%         ax_lim_max = max([abs(ax.XLim), abs(ax.YLim)]);
%         ax.XLim = [-ax_lim_max, ax_lim_max];
%         ax.YLim = [-ax_lim_max, ax_lim_max];
%         hold off;
%         subplot(2,2,4);
%         pbaspect([1 1 1])
%         title('YZ-Plane (Side)');
%         xlabel('Yh');
%         ylabel('Zh');
%         hold on;
%         scatter(mag_calib(2,:), mag_calib(3,:), circle_size, circle_color, 'DisplayName', 'soft/hard iron compensation');
%         plot(axis1(2, :), axis1(3, :), 'r','HandleVisibility','off');
%         plot(axis2(2, :), axis2(3, :), 'g','HandleVisibility','off');
%         plot(axis3(3, :), axis3(3, :), 'b','HandleVisibility','off');
%         % make that all axis have same limit
%         ax = gca; % current axes
%         ax_lim_max = max([abs(ax.XLim), abs(ax.YLim)]);
%         ax.XLim = [-ax_lim_max, ax_lim_max];
%         ax.YLim = [-ax_lim_max, ax_lim_max];
%         hold off;
    
	end
    

end
