% Import simulation Data struct "kiteTruth" which contains:
% - time_s: Timestamp
% - pe_m: inertial position 3D
% - eul_deg: euler angle 3D
% - wb_rps: body angular velocity 3D
% The sensor datas were generated from this data with additional noise and
% bias

% importSimulationData('../../SimulationenBauer/ReferenceData/kite_ref_data_fbauer.mat');

% assumption: only positive time possible
function [kiteTruth, startTime, stopTime] = importSimulationData(file, startTime, stopTime)

load(file);

if ~exist('kiteTruth')
    error('No variable named "kiteTruth"');
end

if startTime < 0
    startTime = kiteTruth.time_s(1);
end

if stopTime < 0
    stopTime = kiteTruth.time_s(end);
end

if startTime > stopTime
    error(['Start time (', startTime, ') cannot be greater than the end time (', stopTime, ')!']);
end

[row_startTime,col,v] = find(kiteTruth.time_s >= startTime,1);
[row_stopTime,col,v] = find(kiteTruth.time_s >= stopTime,1);


scaleTime = 100;
kiteTruth.time_s = kiteTruth.time_s; % .* scaleTime;

%% Movement model
velocity =  [gradient(kiteTruth.pe_m(:,1),kiteTruth.time_s) ...
              gradient(kiteTruth.pe_m(:,2),kiteTruth.time_s) ...
              gradient(kiteTruth.pe_m(:,3),kiteTruth.time_s)]; % [m/s]         
acceleration = [gradient(velocity(:,1),kiteTruth.time_s) ...
              gradient(velocity(:,2),kiteTruth.time_s) ...
              gradient(velocity(:,3),kiteTruth.time_s)]; %[m/s^2]
          
% Assume the first value for accleration is zero to calulate initial values for
% euler angles
acceleration(1,:) = [0 0 0];

vecnorm(acceleration');

% compatible to already existing simulation convention
kiteTruth.time_s = kiteTruth.time_s';
kiteTruth.pe_m = kiteTruth.pe_m';
kiteTruth.eul_deg = kiteTruth.eul_deg';
kiteTruth.wb_rps = kiteTruth.wb_rps';
kiteTruth.ve_mps = velocity';
kiteTruth.ae_mps2 = acceleration';

% % remove first 2 seconds
length_time_s = length(kiteTruth.time_s);
kiteTruth.time_s = (kiteTruth.time_s(:, row_startTime:row_stopTime)-kiteTruth.time_s(row_startTime));
kiteTruth.pe_m = kiteTruth.pe_m(:, row_startTime:row_stopTime);
kiteTruth.eul_deg = kiteTruth.eul_deg(:, row_startTime:row_stopTime);
kiteTruth.wb_rps = kiteTruth.wb_rps(:, row_startTime:row_stopTime);
kiteTruth.ve_mps = kiteTruth.ve_mps(:, row_startTime:row_stopTime);
kiteTruth.ae_mps2 = kiteTruth.ae_mps2(:, row_startTime:row_stopTime);
end