% dbstop if warning
% dbclear if warning
% dbclear all
close all;
try % don't know why it will not closed automatically
    close(figSlider);
catch
end
clc; clear;
close all;

tic;

% load settings of the simulation from different file
AtttitudeEstimationSettings;
                      
 % TODO: 
 % - testing, when gyro offset is calculated at the beginning and then no
 % change
 
%% ##################### Important Notes ##################################
%##########################################################################
% - When the magnetometer is not calibrated, and turning magnetometer off, 
% the roll and pitch estimations are much better. Yaw makes no sense. 
% Finding a method to calibrate magnetometer, use magnetometer only for 
% yaw(so estimating roll and pitch and with the magnetometer the yaw angle 
% instead directly input the magnetometer into the Kalmanfilter), or 
% another solution
% - if the bias is on, make sure P0 of the filters is higher than the bias,
% otherwise the filter is unstable otherwise estimate the bias previously
% and set it as initial value for the filter.
% - determine the inital states carefully otherwise the filter is unstable,
% or set the covariance matrix P0 to really high values, which the initial
% states never exceed!
 
%% Simulation settings
TA = 1; % dummy
if ~importedSystem
    fA = 200; % sampling frequency [Hz]
    TA = 1/fA; % sampling period [s]
    TA_noise = 0; %TA/1000;
    startTime = 0;
    stopTime = 200; % [s]
end

%% Sample rates
if simulationSampleRateForSensors
    TA_sensor.magnetometer = TA/2; % [s]
    TA_sensor.accelerometer = TA/2; % [s]
    TA_sensor.gyro = TA/2; % [s]
    TA_sensor.gps = TA/2; % [s]
    TA_sensor.baro = TA/2; % [s]
    TA_sensor.line = TA/2; % [s]
else
    % sample rates of different sensors
    TA_sensor.magnetometer = 1/200; % [s]
    TA_sensor.accelerometer = 1/200; % [s]
    TA_sensor.gyro = 1/2000; % [s]
    TA_sensor.gps = 1/10; % [s]
    TA_sensor.baro = 1/101; % [s]
    TA_sensor.line = 1/200; % [s]
end

%% Constants
    g_mps2 = 9.81; % gravity [m/s^2]

    % Magnetometer
    if incDec
        % for munich:
        delta_d = 3.12; % [°] declination
        % At +64° the performance is much better
        delta_i = -64.20; % [°] inclination
        if delta_i > 0
            warning('delta_i is wrong! Must be negative');
        end
    else
        warning('Inclination/Declination is set to zero');
        delta_d = 0;
        delta_i = 0;
    end

%% Bias
    if bias_gyro
        warning(['Gyro bias on, so be care that the entries of P0 are high', ...
                'enough to compensate these biases']);
        b_w = [1;2;3]; % bias gyroscope [rad/s]
    else
        b_w = [0;0;0]; % bias gyroscope [rad/s]
    end

    if bias_acc
        b_a = [1;1;1]; % bias accelerometer [m/s^2]
    else
        b_a = [0;0;0]; % bias accelerometer [m/s^2]
    end
    
%% Noise (Variance)
    if noise
        sigma2.w_b = (1.5e-2)^2*[1;1;1]; % [rad/s] sigma2 == sigma^2: attitude process_noise
        sigma2.a_b = (3.5e-1)^2*[1;1;1]*TA; % [m/s^2]variance of the measurement noise of the accelerometer
        sigma2.b_w = (1e-6*[1; 1; 1]); % [rad] variance of process noise of the gyro bias
        sigma2.b_a = [1;1;1]*1e-6;
        sigma2.mag = 3.2400e-06*TA*ones(3,1); % [rad] variance of process noise of the yaw angle from the magnetometer
        sigma2.quaternion = (0.001*ones(4,1)).^2;
        sigma2.baro = 0.6670^2; % [m]
        sigma2.p_gps = (0.6*10)^2 * ones(3, 1); % [m]
        sigma2.v_gps = 1.5^2 * ones(3,1); % [m/s]
        sigma2.length_line = 0.3^2 * TA; %[m]
        sigma2.gamma_line = (3*pi/180).^2; %(atan(10/150)).^2*[1;1]; % in a distance of 150m maximal 10m of error  % ([5;5]*pi/180).^2; %[rad]
        sigma2.p_line = tetherEstimatePositionVariance(sigma2.gamma_line(1), sigma2.length_line);
        
        if varianceEstimation
            addpath('Calibration');
            arduinoPath = '../../Arduino10DOFSensor/DataAndMatlabImport/';
            addpath(arduinoPath);
            foldername = 'VarianceEstimation'; % name of the folder which data should be used
            sensorData = readArduinoData(fullfile(arduinoPath, foldername), 'config.m', 'received.log', g_mps2);
            sigma_2_est = estimateVariance(sensorData);
            sigma2.w_b = sigma_2_est.w_b;
            sigma2.a = sigma_2_est.a_b;
            
            sigma2.w_b = sigma2.w_b;
            sigma2.a_b = sigma2.a_b*TA; % [m/s^2]variance of the measurement noise of the accelerometer
            %sigma2.b_w = ; % [rad] variance of process noise of the gyro bias
            sigma2.mag = sigma2.w_b; % [rad] variance of process noise of the yaw angle from the magnetometer
            clear sigma_2_est;
        end
        
    else
        warning("Noise is off. Please turn on for any statistical filter");
        sigma2.w_b = [0;0;0];
        sigma2.b_w = [0;0;0];
        sigma2.a_b = [0;0;0];
        sigma2.mag = 0*ones(3,1);
        sigma2.quaternion = 0*ones(4,1);
        sigma2.baro = 0; % [m]
        sigma2.p_gps = 0 * ones(3,1); % [m]
        sigma2.v_gps = 0 * ones(3,1); % [m/s]
        sigma2.p_line = 0 * ones(3,1); % [m]
        sigma2.gamma_line = [0;0];
        sigma2.length_line = 0;
        sigma2.b_a = [0;0;0];
    end
    
    sigma2_system = sigma2; % variances to generate measurement data
    if disable_noise
        sigma2_system.w_b = [0;0;0];
        sigma2_system.b_w = [0;0;0];
        sigma2_system.a_b = [0;0;0];
        sigma2_system.mag = 0*ones(3,1);
        sigma2_system.quaternion = 0*ones(4,1);
        sigma2_system.baro = 0; % [m]
        sigma2_system.p_gps = 0; % [m]
        sigma2_system.v_gps = 0; % [m/s]
        sigma2_system.p_line = 0; % [m]
        sigma2_system.gamma_line = [0;0];
        sigma2_system.length_line = 0;
        sigma2_system.b_a = [0;0;0];
    end

%% Magnetometer Calibration
magOffset = [0;0;0];
if calibratingMagnetometer
            addpath('Calibration');
            addpath('Calibration/Source/ellipsoid_fit/ellipsoid_fit');
            addpath('../../Arduino10DOFSensor/DataAndMatlabImport');
            arduinoPath = '../../Arduino10DOFSensor/DataAndMatlabImport/OutsideNewComplete';
            addpath(arduinoPath);
            foldername = 'MagnetometerCalibration'; % name of the folder which data should be used
            sensorData = readArduinoData(fullfile(arduinoPath, foldername), 'config.m', 'received.log', g_mps2);
            useKalmanFilter = 0; % use kalmanfilter to estimate roll and pitch from accelerometer
            plottingMagnetometerCalibration = 1;
            testFolder = 'Rotation3';
            [offsetMag, softIronMatrix]= magnetometerCalibration(sensorData.m_b_unnormalized, ...
                                                                plottingMagnetometerCalibration);
            
            clear arduinoPath foldername sensorData testFolder;
end

if determineMagnetometerAngles
    addpath('Calibration');
    arduinoPath = '../../Arduino10DOFSensor/DataAndMatlabImport/';
    addpath(arduinoPath);
    foldername = 'NotMovingDirectionNorth2'; % name of the folder which data should be used
    sensorData = readArduinoData(fullfile(arduinoPath, foldername), 'config.m', 'received.log', g_mps2);
    [delta_d_new, delta_i_new] = determineInclinationDeclination(sensorData.meas.m_b_unnormalized - magOffset, sensorData.meas.a_b);
    
    if abs(delta_d_new - delta_d) > 5 || abs(delta_i_new - delta_i) > 10
        warning('The estimated magnetic inclination or declination is not how defined in some references found on the internet for munich (delta_d = 2°, delta_i = 64°)')
    end
    delta_d = delta_d_new;
    delta_i = delta_i_new;
    
    clear delta_d_new delta_i_new;
end

%% System initialization
    v0 = [0;0;0];
    p0 = [0;0;0];
if ~importedSystem  
    numberOfRounds = ceil(stopTime/TA);
    if createSystemSource.circle
        gamma_inertia0 = [0; 0; pi/2];
        v0 = [0;0.1571;0];
        p0 = [5;0;0];
    else
        gamma_inertia0 = [0; 0; 0]; % [rad]
        v0 = [0;0;0]; % [m/s]
        p0 = [0;0;0]; % [m]
    end
    system = SystemSimulation(TA, TA_sensor, delta_i, delta_d, g_mps2, sigma2_system);
    truthdata.gamma_inertia = zeros(3, numberOfRounds);
    
    % Init generators used to create the input values
    A_w_phi = 0; % [rad/s]
    A_w_theta = 0; % [rad/s]
    A_w_psi = 0; % [rad/s]
    f_w_phi = 1; % [Hz]
    f_w_theta = 1; % [Hz]
    f_w_psi = 1; % [Hz]
end

%% System Simulation
if importedSystem
    fn = fieldnames(filterSimulation);
    notAllowedSimulations = {'INSEulerEKF', 'INSEulerEKFDST', 'INSMEKF', 'INSMEKFDSTVarianceChange', 'INSMEKFDST_tetherSag'};
    for k=1:numel(fn)
        name = fn{k};
        value = filterSimulation.(fn{k});
        if ~value
            continue;
        end
        for i = 1 : length(notAllowedSimulations)
            if (strcmp(name, notAllowedSimulations(i)))
                error("Not all simulations are allowed when using imported data, because the data don't have tether simulations");
            end
        end
    end
    
    if importSource.useMultikopterData
        path = 'MultikopterData/AttitudeControl_20190727_162819.cdf';
        [measurements, truthdata] = readMultikopterData(path); % truthdata.gamma_inertia = estimated from the copter
        time = measurements.meas.time;
        numberOfRounds = length(time);
        TA = measurements.meas.TA;
        if strcmp(path, 'MultikopterData/AttitudeControl_20190727_162819.cdf')
            % multikopter moved into this position
            % reference position
            truthdata.gamma_inertia = zeros(3, length(time));
        end
        stopTime = time(length(time));
    elseif importSource.useArduinoLog
        arduinoPath = '../../Arduino10DOFSensor/DataAndMatlabImport';
        addpath(arduinoPath);
        foldername = 'OutsideNewComplete/Rotation3';
        [measurements, TA, time_s] = readArduinoData(fullfile(arduinoPath, foldername), 'config.m', 'received.log', g_mps2);
        measurements.meas.m_b = measurements.meas.m_b_unnormalized;
        if exist('offsetMag')
            measurements.meas.m_b_HI_comp = measurements.meas.m_b - offsetMag; % because offset is calculated from this values
            if exist('softIronMatrix')
                %display('Soft Iron calibration is turned off');
                measurements.meas.m_b = softIronMatrix*measurements.meas.m_b_HI_comp;
            end
        end
        m_b_meas_woutOffsetCompensation = measurements.meas.m_b_unnormalized;
        measurements.meas.m_b = measurements.meas.m_b./vecnorm(measurements.meas.m_b);
        time = time_s;
        stopTime = time(length(time));
        numberOfRounds = length(time);
        measurements.meas.p_gps = zeros(3, numberOfRounds);
        measurements.meas.v_gps = zeros(3, numberOfRounds);
        measurements.meas.baro = zeros(1, numberOfRounds);
        measurements.meas.p_line = zeros(3, numberOfRounds);
        truthdata.gamma_inertia = zeros(3, numberOfRounds); % no reference data
        truthdata.v_i = zeros(3, numberOfRounds);
        truthdata.p_i = zeros(3, numberOfRounds);

        clear time_s;
    end
else
    % creating the sensor data from existing inertia acceleration a_i and
    % inertia angular rates w_i
    time = (1:numberOfRounds) * TA;
    ones_time = ones(1, length(time));
    if createSystemSource.circle
        time_to_rotate = 20;
        % radius = 5m
        % fly circle (to test centrifugal acceleration)
        % testing if centrifugal compensation work
        v0 = [0;-((10*pi*cos((2*pi*time_to_rotate/2)/time_to_rotate))/time_to_rotate - (10*pi*cos((2*pi*0)/time_to_rotate))/time_to_rotate)/2;0];
        p0 = [5;0;0];
        gamma_inertia0 = [0;0;-pi/2];
        w_i = [0*time;0*time;2*pi/time_to_rotate*ones_time];
		a_i = [-(20*pi^2*cos((2*pi*time)/time_to_rotate))/time_to_rotate^2;
              -(20*pi^2*sin((2*pi*time)/time_to_rotate))/time_to_rotate^2;
              0*ones_time]; % [m/s^2]
    elseif createSystemSource.stepResponse
        a_i = zeros(3, numberOfRounds);
        timeToReachNewAngle = 0.01; % [s] because if this number is zero, w_i goes to infinity
        % angle at the end is not perfect this due to the rounding below
        phi_old = 0;    % [°]
        phi_new = 10;   % [°]
        theta_old = 0;  % [°]
        theta_new = 10; % [°]
        psi_old = 0;    % [°]
        psi_new = 10;   % [°]
        w_i = zeros(3, numberOfRounds);
        w_i(:, numberOfRounds/2:numberOfRounds/2+round(timeToReachNewAngle/TA)) = deg2rad([phi_new - phi_old; theta_new - theta_old; psi_new - psi_old])/timeToReachNewAngle.*ones(3,round(timeToReachNewAngle/TA)+1);
        gamma_inertia0 = [0;0;0];
        clear timeToReachNewAngle;
    else
        a_i = zeros(3, numberOfRounds);
        w_i = [A_w_phi*sin(2*pi*f_w_phi * time);
               A_w_theta*sin(2*pi*f_w_theta * time);
               A_w_psi*sin(2*pi*f_w_psi * time)];
    end
    disableMagnetometer = 0;
    if disableMagnetometer
        warning('Magnetometer is in the system disabled');
    end
    [measurements, truthdata] = system.simulate(time, a_i, w_i, gamma_inertia0, p0, v0, b_w, disableMagnetometer, sensorFailure);
    gamma_inertia_quaternion = eulerAnglesToQuaternion(truthdata.gamma_inertia);
    truthdata.gamma_inertia = rad2deg(wrapToPi(truthdata.gamma_inertia));
    
    %vecnorm(crossproductMatrix(w_i)*v); %should be constant for a circle; centripetalacceleration
end

if createSliderFigureBool
    SliderFigure; %(time, p, truthdata.gamma_inertia);
end

%% Debug

if debug
    f = figure('Name','Raw Data','NumberTitle','off');
    figure(f);
    movegui(f,'center');
    subplot(4,1,1);
    title('Acceleration in body coords [g]');
    hold('on');
    plot(time, vecnorm(measurements.meas.a_b), 'Displayname', 'Norm');
    plot(time, measurements.meas.a_b(1,:), 'DisplayName', 'x');
    plot(time, measurements.meas.a_b(2,:), 'DisplayName', 'y');
    plot(time, measurements.meas.a_b(3,:), 'DisplayName', 'z');
    hold('off');
    subplot(4,1,2);
    title('Angular rate in body coords [rad/s]');
    hold('on');
    plot(time, measurements.meas.w_b(1,:), 'DisplayName', 'x');
    plot(time, measurements.meas.w_b(2,:), 'DisplayName', 'y');
    plot(time, measurements.meas.w_b(3,:), 'DisplayName', 'z');
    hold('off');
    subplot(4,1,3);
    title('Magnetic field in body coords [-]');
    hold('on');
    plot(time, measurements.meas.m_b(1,:), 'DisplayName', 'x');
    plot(time, measurements.meas.m_b(2,:), 'DisplayName', 'y');
    plot(time, measurements.meas.m_b(3,:), 'DisplayName', 'z');
    hold('off');
    
    a_n = vecnorm(measurements.meas.a_b);
    ax = measurements.meas.a_b(1,:)./a_n;
    ay = measurements.meas.a_b(2,:)./a_n;
    az = measurements.meas.a_b(3,:)./a_n;

    phi_log= atan2(ay,az);
    theta_log= asin(-ax/1);
    % use magnetometer only for yaw angle
    m_n = vecnorm(measurements.meas.m_b);
    mx = measurements.meas.m_b(1,:)./m_n;
    my = measurements.meas.m_b(2,:)./m_n;
    mz = measurements.meas.m_b(3,:)./m_n;
    psi_log =  atan2(sin(phi_log).*mz - cos(phi_log).*my, cos(theta_log).*mx + sin(phi_log).*sin(theta_log).*my + cos(phi_log).*sin(theta_log).*mz);
    subplot(4,1,4);
    title('CalculatedAngles [°]');
    hold('on');
    plot(time, rad2deg(phi_log), 'DisplayName', 'phi');
    plot(time, rad2deg(theta_log), 'DisplayName', 'theta');
    plot(time, rad2deg(psi_log), 'DisplayName', 'psi');
    hold('off');
    
    clear a_n ax ay az phi theta m_n mx my mz psi;
    try
        p_line_calc_sag = zeros(3, length(measurements.meas.length_sag_line));
        for i=1:length(measurements.meas.length_sag_line)
            [x,y,z, ~, ~, ~] = calculatePositionFromLineAngle(measurements.meas.theta_base_line(i), ...
                                                    measurements.meas.theta_kite_line(i), ...
                                                    measurements.meas.phi_sag_line(i), ...
                                                    measurements.meas.length_sag_line(i));
            p_line_calc_sag(:, i) = [x;y;z];
        end
    catch
    end
    
    if exist('truthdata.p_i')
        f = figure('Name', 'Raw Position and velocity');
        figure(f); % must be done, because otherwise the plots are not always in the same figure
        subplot(2,1,1);
        title('Position [m]');
        hold('on');
        if exist('kiteTruth') && isfield(kiteTruth,'pe_m')
            plot(time, truthdata.p_i(1, :), 'DisplayName', 'realX');
            plot(time, truthdata.p_i(2, :), 'DisplayName', 'realY');
            plot(time, truthdata.p_i(3, :), 'DisplayName', 'realZ');
        end
        plot(time, measurements.meas.p_gps(1, :), 'DisplayName', 'gpsX');
        plot(time, measurements.meas.p_gps(2, :), 'DisplayName', 'gpsY');
        plot(time, measurements.meas.p_gps(3, :), 'DisplayName', 'gpsZ');
        plot(time, measurements.meas.baro, 'DisplayName', 'baroZ');
        if exist('p_line_calc_sag') && tetherSag
            plot(time, p_line_calc_sag(1, :), 'c--','DisplayName', 'lineXsag');
            plot(time, p_line_calc_sag(2, :), 'm--','DisplayName', 'lineYsag');
            plot(time, p_line_calc_sag(3, :), 'g--','DisplayName', 'lineZsag');
        end
        plot(time, measurements.meas.p_line(1, :), '-.','DisplayName', 'linePosX');
        plot(time, measurements.meas.p_line(2, :), '-.','DisplayName', 'linePosY');
        plot(time, measurements.meas.p_line(3, :), '-.','DisplayName', 'linePosZ');
        hold('off');
        subplot(2,1,2);
        title('Velocity [m/s]');
        hold('on');
        plot(time, kiteTruth.ve_mps(1, :), 'DisplayName', 'x');
        plot(time, kiteTruth.ve_mps(2, :), 'DisplayName', 'y');
        plot(time, kiteTruth.ve_mps(3, :), 'DisplayName', 'z');
        plot(time, measurements.meas.v_gps(1, :), 'DisplayName', 'gpsX');
        plot(time, measurements.meas.v_gps(2, :), 'DisplayName', 'gpsY');
        plot(time, measurements.meas.v_gps(3, :), 'DisplayName', 'gpsZ');
        hold('off');
    end
    f = figure('Name', 'Sensor availability');
    figure(f);
    title('Show sensor availability');
    minLim = -0.1;
    maxLim = 1.1;
    s(1) = subplot(6, 1, 1);
    stairs(time, measurements.available.acc, 'DisplayName', 'Accelerometer');
    ylim([minLim, maxLim]);
    legend;
    s(2) = subplot(6, 1, 2);
    stairs(time, measurements.available.w_b, 'DisplayName', 'Gyroscope');
    ylim([minLim, maxLim]);
    legend;
    s(3) = subplot(6, 1, 3);
    stairs(time, measurements.available.mag, 'DisplayName', 'Magnetometer');
    ylim([minLim, maxLim]);
    legend;
    s(4) = subplot(6, 1, 4);
    stairs(time, measurements.available.gps, 'DisplayName', 'GPS');
    ylim([minLim, maxLim]);
    legend;
    s(5) = subplot(6, 1, 5);
    stairs(time, measurements.available.baro, 'DisplayName', 'Barometer');
    ylim([minLim, maxLim]);
    legend;
    try
        s(6) = subplot(6, 1, 6);
        stairs(time, measurements.available.line, 'DisplayName', 'Line Angle');
        ylim([minLim, maxLim]);
        legend;
    catch
    end
    
    linkaxes(s, 'x');
    clear p_line_calc s minLim maxLim;
end

%% Filter

if simulation || simulink

%% Filter initialization
    % Initial rotation
    a_n = vecnorm(measurements.meas.a_b);
    a_n = a_n(1);
    % es darf nicht die norm gemacht werden, sonst wird theta immer sehr
    % klein sein!!!
    ax = measurements.meas.a_b(1,1); % ./a_n;
    ay = measurements.meas.a_b(2,1); % ./a_n;
    az = measurements.meas.a_b(3,1); % ./a_n;

    phi_init= atan2(ay,az);
    theta_init= asin(-ax/g_mps2); % hier!
    m_b_init = Rmag_to_e(measurements.meas.m_b(:,1), delta_d, delta_i);
    
    m_n = vecnorm(m_b_init);
    m_n = m_n(1);
    mx = m_b_init(1,1)./m_n;
    my = m_b_init(2,1)./m_n;
    mz = m_b_init(3,1)./m_n;
    psi_init =  atan((sin(phi_init).*mz - cos(phi_init).*my)/ (cos(theta_init).*mx + sin(phi_init).*sin(theta_init).*my + cos(phi_init).*sin(theta_init).*mz));
    
    p_init = measurements.meas.p_gps(:, 1);
    v_init = measurements.meas.v_gps(:, 1);
    
    if filterSimulation.eulerEKF
        x0 = [phi_init; theta_init; psi_init;0;0;0];
        % use high initial covarianzes and the filter converges much faster
        large_angle_uncertainty_rad = 30*pi/180;
        P0 = diag([...
                    [1 1 1]*large_angle_uncertainty_rad^2 ...     % init Euler angle (NED-to-body) uncertainties, rad
                    [1 1 1]*1e-4' ...        % init XYZ gyro bias uncertainties, rad/s
                ]);
        ahrsEuler = AHRSEulerEKF(TA, delta_i, delta_d, g_mps2, x0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
        gamma_i_ahrsEuler = zeros(3, numberOfRounds);
        b_w_ahrsEuler = zeros(3, numberOfRounds);
        
        clear x0 large_angle_uncertainty_rad P0;
    end
    
    if filterSimulation.eulerEKFHIComp
        x0_eulerEKFHIComp = [phi_init; theta_init; psi_init; 0; 0; 0; 0; 0; 0; mx; my; mz;];
        % use high initial covarianzes and the filter converges much faster
        large_angle_uncertainty_rad = 30*pi/180;
        P0_eulerEKFHIComp = diag([...
                    [1 1 1]*large_angle_uncertainty_rad^2 ...     % init Euler angle (NED-to-body) uncertainties, rad
                    [1 1 1]*1e-4 ...        % init XYZ gyro bias uncertainties, rad/s
                    ones(1,3)*1e-1 ...      % init xyz m_e uncertainties
                    ones(1,3)*1e-1          % init xyz m_b uncertainties
                ]);
        ahrsEulerHIComp = AHRSEulerEKFHIComp();
        ahrsEulerHIComp.initFilter(TA, delta_i, delta_d, g_mps2, x0_eulerEKFHIComp, P0_eulerEKFHIComp, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
        gamma_i_ahrsEulerHIComp = zeros(3, numberOfRounds);
        b_w_ahrsEulerHIComp = zeros(3, numberOfRounds);
        
        clear x0 large_angle_uncertainty_rad P0;
    end
    
    if filterSimulation.eulerEKFWoutGyroBias
        x0 = [phi_init; theta_init; psi_init; 0; 0; 0];
        % use high initial covarianzes and the filter converges much faster
        large_angle_uncertainty_rad = 30*pi/180;
        P0 = diag([...
                    [1 1 1]*large_angle_uncertainty_rad^2, ...     % init Euler angle (NED-to-body) uncertainties, rad
                    [1 1 1]
                ]);
        ahrsEulerWoutGyroBias = AHRSEulerEKFWoutGyroBias(TA, delta_i, delta_d, g_mps2, x0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
        gamma_i_ahrsEulerWoutGyroBias = zeros(3, numberOfRounds);
        b_w_ahrsEulerWoutGyroBias = zeros(3, numberOfRounds);
        clear x0 large_angle_uncertainty_rad P0;
    end
   
    if filterSimulation.eulerEKFDifferentSampleTime
        x0 = [phi_init; theta_init; psi_init;0;0;0];
        % use high initial covarianzes and the filter converges much faster
        large_angle_uncertainty_rad = 30*pi/180;
        P0 = diag([...
                [1 1 1]*large_angle_uncertainty_rad^2 ...     % init Euler angle (NED-to-body) uncertainties, rad
                [1 1 1]*1e-4' ...        % init XYZ gyro bias uncertainties, rad/s
                ]);
        ahrsEulerDifferentSampleTime = AHRSEulerEKFDifferentSampleTime(TA, delta_i, delta_d, g_mps2, x0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
        gamma_i_ahrsEulerDifferentSampleTime = zeros(3, numberOfRounds);
        b_w_ahrsEulerDifferentSampleTime = zeros(3, numberOfRounds);
    
        clear x0 large_angle_uncertainty_rad P0;
    end

    if filterSimulation.quaternionEKF
         % norm(q) must be one!
         q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
         x0_quaternionEKF = [q0;0;0;0];
         P0_quaternionEKF = [0.1,0,0,0,0,0,0;
                             0,0.1,0,0,0,0,0;
                             0,0,0.1,0,0,0,0;
                             0,0,0,0.1,0,0,0;
                             0,0,0,0,5,0,0;
                             0,0,0,0,0,5,0;
                             0,0,0,0,0,0,5];
        ahrsQuaternion = AHRSQuaternionEKF();
        ahrsQuaternion.initFilter(TA, x0_quaternionEKF, P0_quaternionEKF, delta_i, delta_d, g_mps2, sigma2.quaternion, sigma2.a_b, sigma2.b_w, sigma2.mag);
        gamma_i_ahrsQuaternion = zeros(3, numberOfRounds);
        b_w_ahrsQuaternion = zeros(3, numberOfRounds);
    
        if (~simulink) % don't delete, because simulink needs this parameter
            clear x0_quaternionEKF P0_quaternionEKF;
        end
    end
    
    if filterSimulation.quaternionEKFDifferentSampleTime
         q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
         x0_quaternionEKF = [q0;0;0;0];
         P0_quaternionEKF = [0.5,0,0,0,0,0,0;
                             0,0.5,0,0,0,0,0;
                             0,0,0.5,0,0,0,0;
                             0,0,0,0.5,0,0,0;
                             0,0,0,0,5,0,0;
                             0,0,0,0,0,5,0;
                             0,0,0,0,0,0,5];
        ahrsQuaternionDifferentSampleTime = AHRSQuaternionEKFDifferentSampleTime();
        ahrsQuaternionDifferentSampleTime.initFilter(TA, x0_quaternionEKF, P0_quaternionEKF, delta_i, delta_d, g_mps2, sigma2.quaternion, sigma2.a_b, sigma2.b_w, sigma2.mag);
        gamma_i_ahrsQuaternionDifferentSampleTime = zeros(3, numberOfRounds);
        b_w_ahrsQuaternionDifferentSampleTime = zeros(3, numberOfRounds);
    
        clear q0 x0_quaternionEKF P0_quaternionEKF;
    end

    if filterSimulation.MEKF
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        P0 = [1, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0;
              0, 0, 0, 5, 0, 0; % bias to a high variance. In real it will never so high. But above I use as bias 3rad/s
              0, 0, 0, 0, 5, 0;
              0, 0, 0, 0, 0, 5];
        gamma_i_ahrsMEKF = zeros(3, numberOfRounds);
        b_w_ahrsMEKF = zeros(3, numberOfRounds); 
        ahrsMEKF = AHRSMEKF();
        ahrsMEKF.initFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
    
        clear q0 b_w0 P0;
    end
    
    if filterSimulation.MEKFDifferentSampleTime
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        gamma_i_ahrsMEKFDifferentSampleTime = zeros(3, numberOfRounds);
        b_w_ahrsMEKFDifferentSampleTime = zeros(3, numberOfRounds); 
        P0 = [1, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0;
              0, 0, 0, 5, 0, 0; % bias to a high variance. In real it will never so high. But above I use as bias 3rad/s
              0, 0, 0, 0, 5, 0;
              0, 0, 0, 0, 0, 5];
        ahrsMEKFDifferentSampleTime = AHRSMEKFDifferentSampleTime(TA, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
    
        clear q0 b_w0 P0;
    end
    
    if filterSimulation.MEKFDifferentSampleTimeVarianceChange
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange = zeros(3, numberOfRounds);
        b_w_ahrsMEKFDifferentSampleTimeVarianceChange = zeros(3, numberOfRounds); 
        P0 = [1, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0;
              0, 0, 0, 5, 0, 0; % bias to a high variance. In real it will never so high. But above I use as bias 3rad/s
              0, 0, 0, 0, 5, 0;
              0, 0, 0, 0, 0, 5];
        ahrsMEKFDifferentSampleTimeVarianceChange = AHRSMEKFDifferentSampleTimeVarianceChange(TA, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
    
        clear q0 b_w0 P0;
    end
    
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        P0 = [1, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0;
              0, 0, 0, 5, 0, 0; % bias to a high variance. In real it will never so high. But above I use as bias 3rad/s
              0, 0, 0, 0, 5, 0;
              0, 0, 0, 0, 0, 5];
        gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer = zeros(3, numberOfRounds);
        b_w_ahrsMEKFDifferentSampleTimeMagnetometer = zeros(3, numberOfRounds); 
        ahrsMEKFDifferentSampleTimeMagnetometer = AHRSMEKFDifferentSampleTimeMagnetometer(TA, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
    
        clear q0 b_w0 P0;
    end
    
    if filterSimulation.MEKFDifferentSampleTimeMagnetometer2
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        P0 = [1, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0;
              0, 0, 0, 5, 0, 0; % bias to a high variance. In real it will never so high. But above I use as bias 3rad/s
              0, 0, 0, 0, 5, 0;
              0, 0, 0, 0, 0, 5];
        gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2 = zeros(3, numberOfRounds);
        b_w_ahrsMEKFDifferentSampleTimeMagnetometer2 = zeros(3, numberOfRounds); 
        ahrsMEKFDifferentSampleTimeMagnetometer2 = AHRSMEKFDifferentSampleTimeMagnetometer2();
        ahrsMEKFDifferentSampleTimeMagnetometer2.initFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
    
        clear q0 b_w0 P0;
    end
    
    if filterSimulation.MEKFwithoutCentripetalCompensation
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        P0 = [1, 0, 0, 0, 0, 0;
              0, 1, 0, 0, 0, 0;
              0, 0, 1, 0, 0, 0;
              0, 0, 0, 5, 0, 0; % bias to a high variance. In real it will never so high. But above I use as bias 3rad/s
              0, 0, 0, 0, 5, 0;
              0, 0, 0, 0, 0, 5];
        gamma_i_ahrsMEKFwithoutCentripetalCompensation = zeros(3, numberOfRounds);
        b_w_ahrsMEKFwithoutCentripetalCompensation = zeros(3, numberOfRounds); 
        ahrsMEKFwithoutCentripetalCompensation = AHRSMEKFwithoutCentripetalCompensation(TA, delta_i, delta_d, g_mps2, q0, b_w0, P0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
    
        clear q0 b_w0 P0;
    end
    
    if filterSimulation.mahonyComplementaryFilter
        % MahonyComplementaryFilter
        q0 = [1;0;0;0];
        b_w0 = [0;0;0];
        mahonyComplementartyFilter = MahonyComplementaryFilter(TA);
        gamma_i_mahonyComplementaryFilter = zeros(3, numberOfRounds);
        
        clear q0 b_w0;
    end
    
    if filterSimulation.mahonyExplicitComplementaryFilter
        % MahonyExplicitComplementaryFilter
        q0 = [1;0;0;0];
        b_w0 = [0;0;0];
        mahonyExplicitComplementaryFilter = MahonyExplicitComplementaryFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0);
        
        gamma_i_mahonyExplicitComplementaryFilter = zeros(3, numberOfRounds);
        b_w_mahonyExplicitComplementaryFilter = zeros(3, numberOfRounds);
        
        clear q0 b_w0;
    end
    
	if filterSimulation.INSEulerEKF
        % INS EulerEKF
        x0 = [phi_init;theta_init;psi_init;p_init;v_init;0;0;0];
        %x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
        % use high initial covarianzes and the filter converges much faster
        large_angle_uncertainty_rad = 30*pi/180;
        P0 = diag([...
        [1 1 1]*large_angle_uncertainty_rad^2 ...     % init Euler angle (NED-to-body) uncertainties, rad
        [1 1 1]*10 ... % position
        [1 1 1]*10 ... % velocity
        [1 1 1]*5 ...        % init XYZ gyro bias uncertainties, rad/s
        ]);
        %sigma2_process_noise = diag([sigma2.w_b', sigma2.v_gps', sigma2.a_b', sigma2.b_w']);
        %sigma2_measurement_noise = diag([sigma2.mag', sigma2.baro, sigma2.p_gps', sigma2.p_line', sigma2.v_gps']);
        sigma2_process_noise = diag([sigma2.w_b', sigma2.v_gps', sigma2.a_b', sigma2.b_w']);
        sigma2_measurement_noise = diag([sigma2.mag', sigma2.baro, sigma2.p_gps', sigma2.p_line', sigma2.v_gps']);
        insEulerEKF = INSEulerEKF(TA, delta_i, delta_d, g_mps2, x0, P0, sigma2_process_noise, sigma2_measurement_noise);
        gamma_i_insEulerEKF = zeros(3, numberOfRounds);
        b_w_insEulerEKF = zeros(3, numberOfRounds);
        p_insEulerEKF = zeros(3, numberOfRounds);
        v_insEulerEKF = zeros(3, numberOfRounds);
        
        clear sigma2_process_noise sigma2_measurement_noise P0 x0;
    end
    
	if filterSimulation.INSEulerEKFDST
        % INS EulerEKF
        x0 = [phi_init;theta_init;psi_init;p_init;v_init;0;0;0];
        % use high initial covarianzes and the filter converges much faster
        large_angle_uncertainty_rad = 30*pi/180;
        P0 = diag([...
        [1 1 1]*large_angle_uncertainty_rad^2 ...     % init Euler angle (NED-to-body) uncertainties, rad
        [1 1 1]*10 ... % position
        [1 1 1]*5 ... % velocity
        [1 1 1]*10 ...        % init XYZ gyro bias uncertainties, rad/s
        ]);
        sigma2_process_noise = diag([sigma2.w_b', sigma2.v_gps', sigma2.a_b', sigma2.b_w']);
        sigma2_measurement_noise = diag([sigma2.mag', sigma2.baro, sigma2.p_gps', sigma2.p_line', sigma2.v_gps']);
        insEulerEKFDST = INSEulerEKFDST(TA, delta_i, delta_d, g_mps2, x0, P0, sigma2_process_noise, sigma2_measurement_noise);
        gamma_i_insEulerEKFDST = zeros(3, numberOfRounds);
        b_w_insEulerEKFDST = zeros(3, numberOfRounds);
        p_insEulerEKFDST = zeros(3, numberOfRounds);
        v_insEulerEKFDST = zeros(3, numberOfRounds);
        clear sigma2_process_noise sigma2_measurement_noise P0 x0;
    end
    
    if filterSimulation.INSMEKF
         % AHRS QuaternionenMEKFw/outCC
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        x0 = [0;0;0;p_init;v_init;0;0;0];
        P0 = diag([...
        [1 1 1]*1 ...    % angle error (a)
        [1 1 1]*3 ... % position
        [1 1 1]*3 ... % velocity
        [1 1 1]*100 ...        % init XYZ gyro bias uncertainties, rad/s
        ]);
        gamma_i_insMEKF = zeros(3, numberOfRounds);
        p_insMEKF = zeros(3, numberOfRounds);
        v_insMEKF = zeros(3, numberOfRounds);
        b_w_insMEKF = zeros(3, numberOfRounds); 
        insMEKF = INSMEKF();
        sigma2_process_noise = diag([sigma2.w_b', sigma2.v_gps', sigma2.a_b', sigma2.b_w']);
        sigma2_measurement_noise = diag([sigma2.mag', sigma2.baro, sigma2.p_gps', sigma2.p_line', sigma2.v_gps']);
        insMEKF.initFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0, x0, P0, sigma2_process_noise, sigma2_measurement_noise);
    
        clear q0 b_w0 sigma2_process_noise sigma2_measurement_noise P0 x0;
    end
    
     if filterSimulation.INSMEKFDST
         % AHRS QuaternionenMEKFw/outCC
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        x0 = [0;0;0;p_init;v_init;b_w0];
        P0 = diag([...
        [1 1 1]*1 ...    % angle error (a)
        [1 1 1]*3 ... % position
        [1 1 1]*3 ... % velocity
        [1 1 1]*10^2 ...        % init XYZ gyro bias uncertainties, rad/s
        ]);
        gamma_i_INSMEKFDST = zeros(3, numberOfRounds);
        p_INSMEKFDST = zeros(3, numberOfRounds);
        v_INSMEKFDST = zeros(3, numberOfRounds);
        b_w_INSMEKFDST = zeros(3, numberOfRounds); 
        INSMEKFDST = INSMEKFDST();
        sigma2_process_noise = [sigma2.w_b; sigma2.v_gps; sigma2.a_b; sigma2.b_w];
        sigma2_measurement_noise = [sigma2.mag; sigma2.baro; sigma2.p_gps; sigma2.p_line; sigma2.v_gps];
        INSMEKFDST.initFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0, x0, P0, sigma2_process_noise, sigma2_measurement_noise);
    
        clear q0 b_w0 sigma2_process_noise sigma2_measurement_noise P0 x0;
     end
    
    if filterSimulation.INSMEKFDSTVarianceChange
         % AHRS QuaternionenMEKFw/outCC
        q0 = eulerAnglesToQuaternion([phi_init; theta_init; psi_init]);
        b_w0 = [0;0;0];
        x0 = [0;0;0;p_init;v_init;0;0;0];
        P0 = diag([...
        [1 1 1]*1 ...    % angle error (a)
        [1 1 1]*3 ... % position
        [1 1 1]*3 ... % velocity
        [1 1 1]*100 ...        % init XYZ gyro bias uncertainties, rad/s
        ]);
        gamma_i_INSMEKFDSTVarianceChange = zeros(3, numberOfRounds);
        p_INSMEKFDSTVarianceChange = zeros(3, numberOfRounds);
        v_INSMEKFDSTVarianceChange = zeros(3, numberOfRounds);
        b_w_INSMEKFDSTVarianceChange = zeros(3, numberOfRounds); 
        INSMEKFDSTVarianceChange = INSMEKFDSTVarianceChange();
        sigma2_process_noise = [sigma2.w_b; sigma2.v_gps; sigma2.a_b; sigma2.b_w];
        sigma2_measurement_noise = [sigma2.mag; sigma2.baro; sigma2.p_gps; sigma2.p_line; sigma2.v_gps];
        INSMEKFDSTVarianceChange.initFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0, x0, P0, sigma2_process_noise, sigma2_measurement_noise);
    
        clear q0 b_w0 sigma2_process_noise sigma2_measurement_noise P0 x0;
    end
    
    if filterSimulation.INSMEKFDST_tetherSag
         % AHRS QuaternionenMEKFw/outCC
        q0 = eulerAnglesToQuaternion([0; theta_init; psi_init]);
        b_w0 = [0;0;0];
        x0 = [0;0;0;p_init;v_init;0;0;0];
        P0 = diag([...
        [1 1 1]*0.3 ...    % angle error (a)
        [1 1 1]*10 ... % position
        [1 1 1]*5 ... % velocity
        [1 1 1]*10 ...        % init XYZ gyro bias uncertainties, rad/s
        ]);
        gamma_i_INSMEKFDST_tetherSag = zeros(3, numberOfRounds);
        p_INSMEKFDST_tetherSag = zeros(3, numberOfRounds);
        v_INSMEKFDST_tetherSag = zeros(3, numberOfRounds);
        b_w_INSMEKFDST_tetherSag = zeros(3, numberOfRounds); 
        INSMEKFDST_tetherSag = INSMEKFDST_tetherSag();
        sigma2_process_noise = [sigma2.w_b; sigma2.v_gps; sigma2.a_b; sigma2.b_w];
        sigma2_measurement_noise = [sigma2.mag; sigma2.baro; sigma2.p_gps; sigma2.p_line; sigma2.v_gps];
        INSMEKFDST_tetherSag.initFilter(TA, delta_i, delta_d, g_mps2, q0, b_w0, x0, P0, sigma2_process_noise, sigma2_measurement_noise);
    
        clear q0 b_w0 sigma2_process_noise sigma2_measurement_noise P0 x0;
    end
    
    clear a_n ax ay az m_n mx my mz phi_init theta_init psi_init m_b_init;
end

if simulation
%% Filter simulation

    display('Start simulation');
    n=0;
    tic;
    for i=1:numberOfRounds

    %     fprintf(repmat('\b',1,n)); % delete previous message
    %     msg = [num2str(i/numberOfRounds*100), '%%'];
    %     fprintf(msg);
    %     n=numel(msg);

        if filterSimulation.eulerEKF
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsEuler.performEKF(u, y);
            [gamma_i, b_w_hat] = ahrsEuler.states();
            gamma_i_ahrsEuler(:,i) = rad2deg(gamma_i);
            b_w_ahrsEuler(:,i) = b_w_hat;
        end
        
        if filterSimulation.eulerEKFHIComp
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsEulerHIComp.performEKF(u, y);
            [gamma_i, b_w_hat] = ahrsEulerHIComp.states();
            gamma_i_ahrsEulerHIComp(:,i) = rad2deg(gamma_i);
            b_w_ahrsEulerHIComp(:,i) = b_w_hat;
        end
        
        if filterSimulation.eulerEKFWoutGyroBias
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsEulerWoutGyroBias.performEKF(u, y);
            gamma_i = ahrsEulerWoutGyroBias.states();
            gamma_i_ahrsEulerWoutGyroBias(:,i) = rad2deg(gamma_i);
            b_w_ahrsEulerWoutGyroBias(:,i) = ahrsEulerWoutGyroBias.gyroBias();
        end
        
        if filterSimulation.eulerEKFDifferentSampleTime
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsEulerDifferentSampleTime.performEKF(u, y, measurements.available.acc(i), measurements.available.mag(i));
            [gamma_i, b_w_hat] = ahrsEulerDifferentSampleTime.states();
            gamma_i_ahrsEulerDifferentSampleTime(:,i) = rad2deg(gamma_i);
            b_w_ahrsEulerDifferentSampleTime(:,i) = b_w_hat;
        end

        if filterSimulation.quaternionEKF
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsQuaternion.performEKF(u, y);
            [q_hat, b_w_hat] = ahrsQuaternion.states();
            q = Quaternion(q_hat);
            gamma_i_ahrsQuaternion(:,i) = rad2deg(q.quatToEuler());
            %gamma_i_ahrsQuaternionDebug(:,i) = q.debugAngleDeg();
            b_w_ahrsQuaternion(:,i) = b_w_hat;
        end

        if filterSimulation.quaternionEKFDifferentSampleTime
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsQuaternionDifferentSampleTime.performEKF(u, y, measurements.available.acc(i), measurements.available.mag(i));
            [q_hat, b_w_hat] = ahrsQuaternionDifferentSampleTime.states();
            q = Quaternion(q_hat);
            gamma_i_ahrsQuaternionDifferentSampleTime(:,i) = rad2deg(q.quatToEuler());
            b_w_ahrsQuaternionDifferentSampleTime(:,i) = b_w_hat;
        end
        
        if filterSimulation.MEKF
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            K(:,:,i) = ahrsMEKF.performEKF(u, y);
            q_hat = ahrsMEKF.attitude();
            q = Quaternion(q_hat);
            b_w_hat = ahrsMEKF.gyroBias();
            gamma_i_ahrsMEKF(:, i) = rad2deg(q.quatToEuler());
            b_w_ahrsMEKF(:, i) = b_w_hat;
        end
        
        if filterSimulation.MEKFDifferentSampleTime
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsMEKFDifferentSampleTime.performEKF(u, y, measurements.available.acc(i), measurements.available.mag(i));
            q_hat = ahrsMEKFDifferentSampleTime.attitude();
            q = Quaternion(q_hat);
            b_w_hat = ahrsMEKFDifferentSampleTime.gyroBias();
            gamma_i_ahrsMEKFDifferentSampleTime(:, i) = rad2deg(q.quatToEuler());
            b_w_ahrsMEKFDifferentSampleTime(:, i) = b_w_hat;
        end
        
        if filterSimulation.MEKFDifferentSampleTimeVarianceChange
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsMEKFDifferentSampleTimeVarianceChange.performEKF(u, y, measurements.available.acc(i), measurements.available.mag(i));
            q_hat = ahrsMEKFDifferentSampleTimeVarianceChange.attitude();
            q = Quaternion(q_hat);
            b_w_hat = ahrsMEKFDifferentSampleTimeVarianceChange.gyroBias();
            gamma_i_ahrsMEKFDifferentSampleTimeVarianceChange(:, i) = rad2deg(q.quatToEuler());
            b_w_ahrsMEKFDifferentSampleTimeVarianceChange(:, i) = b_w_hat;
        end
        
        if filterSimulation.MEKFDifferentSampleTimeMagnetometer
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsMEKFDifferentSampleTimeMagnetometer.performEKF(u, y, measurements.available.acc(i), measurements.available.mag(i));
            q_hat = ahrsMEKFDifferentSampleTimeMagnetometer.attitude();
            q = Quaternion(q_hat);
            b_w_hat = ahrsMEKFDifferentSampleTimeMagnetometer.gyroBias();
            gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer(:, i) = rad2deg(q.quatToEuler());
            b_w_ahrsMEKFDifferentSampleTimeMagnetometer(:, i) = b_w_hat;
        end
        
        if filterSimulation.MEKFDifferentSampleTimeMagnetometer2
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            ahrsMEKFDifferentSampleTimeMagnetometer2.performEKF(u, y, measurements.available.acc(i), measurements.available.mag(i));
            q_hat = ahrsMEKFDifferentSampleTimeMagnetometer2.attitude();
            q = Quaternion(q_hat);
            b_w_hat = ahrsMEKFDifferentSampleTimeMagnetometer2.gyroBias();
            gamma_i_ahrsMEKFDifferentSampleTimeMagnetometer2(:, i) = rad2deg(q.quatToEuler());
            b_w_ahrsMEKFDifferentSampleTimeMagnetometer2(:, i) = b_w_hat;
        end
        
        if filterSimulation.MEKFwithoutCentripetalCompensation
            u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
            y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
            K(:,:,i) = ahrsMEKFwithoutCentripetalCompensation.performEKF(u, y);
            q_hat = ahrsMEKFwithoutCentripetalCompensation.attitude();
            q = Quaternion(q_hat);
            b_w_hat = ahrsMEKFwithoutCentripetalCompensation.gyroBias();
            gamma_i_ahrsMEKFwithoutCentripetalCompensation(:, i) = rad2deg(q.quatToEuler());
            b_w_ahrsMEKFwithoutCentripetalCompensation(:, i) = b_w_hat;
        end
          
        if filterSimulation.mahonyComplementaryFilter
            mahonyComplementartyFilter.update(measurements.meas.a_b(:,i), measurements.meas.m_b(:,i), measurements.meas.w_b(:,i));
            gamma_i_mahonyComplementaryFilter(:, i) = rad2deg(mahonyComplementartyFilter.eulerAngles());
        end
        
        if filterSimulation.mahonyExplicitComplementaryFilter
            mahonyExplicitComplementaryFilter.update(measurements.meas.a_b(:,i), measurements.meas.m_b(:,i), measurements.meas.w_b(:,i), measurements.meas.p_gps(:, i), measurements.meas.v_gps(:, i));
            gamma_i_mahonyExplicitComplementaryFilter(:, i) = rad2deg(mahonyExplicitComplementaryFilter.eulerAngles());
            b_w_mahonyExplicitComplementaryFilter(:, i) = mahonyExplicitComplementaryFilter.bias();
        end
        
        if filterSimulation.INSEulerEKF
            u = [measurements.meas.w_b(:,i); measurements.meas.a_b(:,i)];
            y = [measurements.meas.m_b(:,i); measurements.meas.baro(i); measurements.meas.p_gps(:,i); measurements.meas.p_line(:,i); measurements.meas.v_gps(:,i)];
            insEulerEKF.performEKF(u, y);
            [gamma_i, pos, vel, b_w_hat] = insEulerEKF.states();
            gamma_i_insEulerEKF(:,i) = rad2deg(gamma_i);
            b_w_insEulerEKF(:,i) = b_w_hat;
            p_insEulerEKF(:, i) = pos;
            v_insEulerEKF(:, i) = vel;
        end
        
        if filterSimulation.INSEulerEKFDST
            u = [measurements.meas.w_b(:,i); measurements.meas.a_b(:,i)];
            y = [measurements.meas.m_b(:,i); measurements.meas.baro(i); measurements.meas.p_gps(:,i); measurements.meas.p_line(:,i); measurements.meas.v_gps(:,i)];
            available = [measurements.available.mag(i)*ones(3,1); measurements.available.baro(i); measurements.available.gps(i)*ones(3,1); measurements.available.line(i)*ones(3,1); measurements.available.gps(i)*ones(3,1)];
            insEulerEKFDST.performEKF(u, y, available);
            [gamma_i, pos, vel, b_w_hat] = insEulerEKFDST.states();
            gamma_i_insEulerEKFDST(:,i) = rad2deg(gamma_i);
            b_w_insEulerEKFDST(:,i) = b_w_hat;
            p_insEulerEKFDST(:, i) = pos;
            v_insEulerEKFDST(:, i) = vel;
        end
        
        if filterSimulation.INSMEKF
            u = [measurements.meas.w_b(:,i); measurements.meas.a_b(:,i)];
            y = [measurements.meas.m_b(:,i); measurements.meas.baro(i); measurements.meas.p_gps(:,i); measurements.meas.p_line(:,i); measurements.meas.v_gps(:,i)];
            insMEKF.performEKF(u, y);
            [gamma_i, pos, vel, b_w_hat] = insMEKF.insStates();
            gamma_i_insMEKF(:,i) = rad2deg(gamma_i);
            b_w_insMEKF(:,i) = b_w_hat;
            p_insMEKF(:, i) = pos;
            v_insMEKF(:, i) = vel;
        end
        
        if filterSimulation.INSMEKFDSTVarianceChange
            u = [measurements.meas.w_b(:,i); measurements.meas.a_b(:,i)];
            y = [measurements.meas.m_b(:,i); measurements.meas.baro(i); measurements.meas.p_gps(:,i); measurements.meas.p_line(:,i); measurements.meas.v_gps(:,i)];
            INSMEKFDSTVarianceChange.performEKF(u, y, measurements.available.acc(i), measurements.available.gyr(i), ... 
                                                      measurements.available.mag(i), measurements.available.baro(i), ...
                                                      measurements.available.gps(i), measurements.available.line(i));
            [gamma_i, pos, vel, b_w_hat] = INSMEKFDSTVarianceChange.insStates();
            gamma_i_INSMEKFDSTVarianceChange(:,i) = rad2deg(gamma_i);
            b_w_INSMEKFDSTVarianceChange(:,i) = b_w_hat;
            p_INSMEKFDSTVarianceChange(:, i) = pos;
            v_INSMEKFDSTVarianceChange(:, i) = vel;
        end
        
        if filterSimulation.INSMEKFDST
            u = [measurements.meas.w_b(:,i); measurements.meas.a_b(:,i)];
            y = [measurements.meas.m_b(:,i); measurements.meas.baro(i); measurements.meas.p_gps(:,i); measurements.meas.p_line(:,i); measurements.meas.v_gps(:,i)];
            INSMEKFDST.performEKF(u, y, measurements.available.acc(i), measurements.available.gyr(i), ... 
                                      measurements.available.mag(i), measurements.available.baro(i), ...
                                      measurements.available.gps(i), measurements.available.line(i));
            [gamma_i, pos, vel, b_w_hat] = INSMEKFDST.insStates();
            gamma_i_INSMEKFDST(:,i) = rad2deg(gamma_i);
            b_w_INSMEKFDST(:,i) = b_w_hat;
            p_INSMEKFDST(:, i) = pos;
            v_INSMEKFDST(:, i) = vel;
        end
        
        if filterSimulation.INSMEKFDST_tetherSag
            u = [measurements.meas.w_b(:,i); measurements.meas.a_b(:,i)];
            line = [measurements.meas.theta_base_line(i);
                    measurements.meas.theta_kite_line(i)
                    measurements.meas.phi_sag_line(i)
                    measurements.meas.length_sag_line(i)];
            
            y = [measurements.meas.m_b(:,i); measurements.meas.baro(i); measurements.meas.p_gps(:,i); line; measurements.meas.v_gps(:,i)];
            INSMEKFDST_tetherSag.performEKF(u, y, measurements.available.acc(i), measurements.available.gyr, ... 
                                                      measurements.available.mag(i), measurements.available.baro(i), ...
                                                      measurements.available.gps(i), measurements.available.line(i));
            [gamma_i, pos, vel, b_w_hat] = INSMEKFDST_tetherSag.insStates();
            gamma_i_INSMEKFDST_tetherSag(:,i) = rad2deg(gamma_i);
            b_w_INSMEKFDST_tetherSag(:,i) = b_w_hat;
            p_INSMEKFDST_tetherSag(:, i) = pos;
            v_INSMEKFDST_tetherSag(:, i) = vel;
        end
        
    end
    clear u y available gamma_i pos vel;
else
    disp('No simulation started (simulation=0)');
end
duration = toc;
disp(['Duration to complete simulation of the filters: ', num2str(duration), 's']);

%% ahrsQuaternionenMEKFDifferentSampleTimeSpeedOptimized
% evaluation: speed improvement is not so high, but maybe in the generated
% code its better
if evalSpeedOptimizedFilter
    warning('Speedoptimized filter will be evaluated');
    if filterSimulation.MEKFDifferentSampleTime
        % AHRS QuaternionenMEKFDST
        q0 = [1;0;0;0];
        b_w0 = [0;0;0];
        gamma_i_ahrsMEKFDifferentSampleTimeSO = zeros(3, numberOfRounds);
        b_w_ahrsMEKFDifferentSampleTimeSO = zeros(3, numberOfRounds); 
        ahrsMEKFDifferentSampleTimeSO = AHRSMEKFDifferentSampleTimeSpeedOptimized(TA, delta_i, delta_d, g_mps2, q0, b_w0, sigma2.w_b, sigma2.a_b, sigma2.b_w, sigma2.mag);
    end
    tic;
    for i=1:numberOfRounds
        % ahrsQuaternionenMEKFDifferentSampleTimeSpeedOptimized
        u = [measurements.meas.w_b(:,i); measurements.meas.v_gps(:,i)];
        y = [measurements.meas.a_b(:,i); measurements.meas.m_b(:,i)];
        ahrsMEKFDifferentSampleTimeSO.performEKF(u, y, measurements.available.acc(i), measurements.available.mag(i));
        q_hat = ahrsMEKFDifferentSampleTimeSO.attitude();
        q = Quaternion(q_hat);
        b_w_hat = ahrsMEKFDifferentSampleTimeSO.gyroBias();
        gamma_i_ahrsMEKFDifferentSampleTimeSO(:, i) = rad2deg(q.quatToEuler());
    end
    optimizedFilter = toc;
    disp(['Optimized filter time: ', num2str(optimizedFilter), 's']);
end

if filterSimulation.MEKF && debug && simulation
    
    display('Plot Kalman Gain of MEKF');
    f = figure('Name','Plot Kalman Gain of MEKF');
    figure(f);
    K_strich = permute(K, [3,2,1]);
    hold('on');
    plot(time, K_strich(:,1,1)');
    plot(time, K_strich(:,1,2)');
    plot(time, K_strich(:,1,3)');
    plot(time, K_strich(:,1,4)');
    plot(time, K_strich(:,1,5)');
    plot(time, K_strich(:,1,6)');
    plot(time, K_strich(:,2,1)');
    plot(time, K_strich(:,2,2)');
    plot(time, K_strich(:,2,3)');
    plot(time, K_strich(:,2,4)');
    plot(time, K_strich(:,2,5)');
    plot(time, K_strich(:,2,6)');
    plot(time, K_strich(:,3,1)');
    plot(time, K_strich(:,3,2)');
    plot(time, K_strich(:,3,3)');
    plot(time, K_strich(:,3,4)');
    plot(time, K_strich(:,3,5)');
    plot(time, K_strich(:,3,6)');
    plot(time, K_strich(:,4,1)');
    plot(time, K_strich(:,4,2)');
    plot(time, K_strich(:,4,3)');
    plot(time, K_strich(:,4,4)');
    plot(time, K_strich(:,4,5)');
    plot(time, K_strich(:,4,6)');
    plot(time, K_strich(:,5,1)');
    plot(time, K_strich(:,5,2)');
    plot(time, K_strich(:,5,3)');
    plot(time, K_strich(:,5,4)');
    plot(time, K_strich(:,5,5)');
    plot(time, K_strich(:,5,6)');
    plot(time, K_strich(:,6,1)');
    plot(time, K_strich(:,6,2)');
    plot(time, K_strich(:,6,3)');
    plot(time, K_strich(:,6,4)');
    plot(time, K_strich(:,6,5)');
    plot(time, K_strich(:,6,6)');
    hold('off')
end

%%
Evaluation;

if simulink || createCCode
    simulink_u = [time', measurements.meas.w_b', measurements.meas.v_gps'];
    simulink_y = [time', measurements.meas.a_b', measurements.meas.m_b'];
    simulink_acc_available = [time', measurements.available.acc'];
    simulink_mag_available = [time', measurements.available.mag'];
end

if createCCode
    file = 'AttitudeEstimation'; % Simulink model
    system_to_compile = 'AHRSQuaternionEKF';
    system_name = struct();
    system_name.name = file;
    system_name.model = system_to_compile;
    
    load_system(system_name.name)
    open_system([system_name.name, '/',system_name.model])
    rtwbuild([system_name.name, '/',system_name.model])
end

if optimize
    display('Start optimization');
    OptimizeFilters;
else
    display('No optimizations (optimize=0)');
end
