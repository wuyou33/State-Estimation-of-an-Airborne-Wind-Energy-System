% this file includes all settings used in main.m
addpath('helpfunctions');
addpath('Filters');
addpath('System');
addpath('DataGeneration'); % read truth data from file and generate Sensor data
addpath('OptimizationFunctions');
addpath('FigureScripts');
addpath('MultikopterData');
addpath('TetherSimulation');
addpath('LatexExport');

%% sources:
% [1]: Barton2012 - Fundamentals of Small Unmanned Aircraft Flight
% [2]: Maley2013  - Multiplicative Quaternion Extended Kalman Filtering for Nonspinning Guided Projectiles
% [3]: Mattia Giurato - https://github.com/Murmele/Attitude_Estimation

plotting = 1; % plotting of different plots. See Evaluation.m
dispAvailable = 1; % 1: display availability of the sensors in plots next to the estimated attitudes
evaluation = 1; % table with the RMSD results
maxErrorTime = 50; % [s] time from which on the maximum error should be found
optimize = 0; % optimize parameters of the filter
simulation = 1; % 0: skip simulation, 1: execute simulation
simulink = 1; % initialize simulink from workspace blocks
createCCode = 0; % create CCode of the Simulink model
evalSpeedOptimizedFilter = 0; % evaluate speed optimized MEKF

% values are quite small. Rechecking function
varianceEstimation = 0; % 1: import data and estimate variance of them; 0: use predefined values
calibratingMagnetometer = 0;
    plotDifferenceAngles = 0; % plots angle calculations, without compensation, with compensation of hard iron effects, with hard and soft iron compensation
determineMagnetometerAngles = 0;

debug = 1; % plots and scripts for debugging executed
%removeTimePlotEvaluation = 1; % Time which should be ignored by the plot and also the evaluation (so that filter can reach steady state)
createSliderFigureBool = 0; % 1: show figure with truth position and truth orientation

startTime = 0; % startTime from the importedData (negative means first value used)
stopTime = -1; % stopTime from the importedData (negative means last value used)
importedSystem = 1; % 0: generate Sensor data from own system. 1: generate Sensor data from imported system
    importSource = struct('useMultikopterData' , 1, ... %1: use data exported from my multikopter. Previouse source must be zero to use this source!
                          'useArduinoLog', 1);          % 1: use exported data from the arduino with 10DOF sensor.  
    % if nothing selected. See in system initialization the behavior
    createSystemSource = struct('circle', 0, ...        % fly a circle to test centrifugal compensation
                                'stepResponse', 0);     % make a angle step to test the performance. circle must be zero!
% in some cases the "DifferentSampleTime" filters are much better
simulationSampleRateForSensors = 0; % 0: use defined sample rates for the sensors. 1: use TA/2 as sample rate for the sensor
bias_gyro = 1; % 1: turn on gyro bias defined below, 0: gyro bias = 0
bias_acc = 0; % not compensated in the most filters
noise = 1; % 1: add noise (must be on for all statistical filters (EKF, MEKF), because otherwise the variances are zero, 0: no noise
disable_noise = 0; % 1: noise in the generated data is zero. 0: noise is added with the variances defined below. It do not disable the covariances, only the noise on the generated signals (works only for importedSystem=0)
incDec = 1; % inclination/ declination of magnetic field. 0: inclination = declination = 0°, 1: see values where this variable is in a condition
disableGPSatToHighAccelerations = 0;
tetherSag = 0; % 1: the line length is about 3% higher than the distance to simulate tether sag. If filter don't simulate the tether sag, there is a offset on the position and the filter is anymore good
loadSimulationData = 0; % if simulation was already one time done, a mat file with the simulation data was created. So for the next time this simulation data is taken
% Select filter which should be simulated:
filterSimulation = struct('eulerEKF', 1, ...% ignores sensor availability
                          'eulerEKFHIComp', 1, ... ignores sensor availability
                          'eulerEKFWoutGyroBias', 0, ...% ignores sensor availability
                          'eulerEKFDifferentSampleTime', 0, ...
                          'quaternionEKF', 1, ... % ignores sensor availability
                          'quaternionEKFDifferentSampleTime', 0, ...
                          'MEKF', 1, ...% ignores sensor availability
                          'MEKFDifferentSampleTime', 0, ...
                          'MEKFDifferentSampleTimeVarianceChange', 0, ...
                          'MEKFDifferentSampleTimeMagnetometer', 0, ...
                          'MEKFDifferentSampleTimeMagnetometer2', 0, ...
                          'MEKFwithoutCentripetalCompensation', 0, ...
                          'mahonyComplementaryFilter',0, ... % not working properly ( no centrifugal compensation)
                          'mahonyExplicitComplementaryFilter', 0, ...% ignores sensor availability
                          'INSEulerEKF', 0, ... % ignores sensor availability
                          'INSEulerEKFDST', 0, ...
                          'INSMEKF', 0, ...
                          'INSMEKFDST', 0, ...
                          'INSMEKFDSTVarianceChange', 0, ...
                          'INSMEKFDST_tetherSag', 0 ...
                          );
angleErrorPlot = 'INSEulerEKF'; % set the filter for which the error plot should be plotted.
gyroBiasPlot = 'INSEulerEKF'; % set the plot data for the gyro bias plot
posVelPlot = 'INSEulerEKF';

% Select filter which should be optimized:
filterOptimization = struct('eulerEKF', 0, ...
                          'quaternionEKF', 0, ...
                          'MEKF',0, ...
                          'mahonyComplementaryFilter',0, ...
                          'mahonyExplicitComplementaryFilter', 1, ...
                          'MEKFDifferentSampleTimeVarianceChange', 0 ...
                          );
 
% Simulating sensor failure
% fail: if true, the sensor will fail at time "time"
% time [s]: time at which the sensor fails and any other values are available
% negative time value means the sensor is anytime available
% setting time to high, the sensor will never fail
sensorFailure.magnetometer.fail = false;
sensorFailure.magnetometer.time = 90;
sensorFailure.accelerometer.fail = false;
sensorFailure.accelerometer.time = 10;
sensorFailure.gyroscope.fail = false;
sensorFailure.gyroscope.time = 100;
sensorFailure.gps.fail = false;
sensorFailure.gps.time = 90;
sensorFailure.lineAngle.fail = false;
sensorFailure.lineAngle.time = 40;
sensorFailure.barometer.fail = false;
sensorFailure.barometer.time = 50;
% if this variable is true, the kite don't move again and just tries to
% land without damage. If this variable is false, the kite is moving his
% trajectory (seeing, if the estimation is quite good)
sensorFailure.landing = true;
                      
createDifferentSampleTimeTable = 0; % creates latex table for different Sample Time
create4gLimitTable = 0; % createx Latex table for 4g GPS limit
createTetherSagTable = 1; % create Latex table for the simulation of tether sag