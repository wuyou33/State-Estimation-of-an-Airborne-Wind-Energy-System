% Function to read data from a csv file (not used, because there are too
% less datas in these csv files provided from florian. Only one complete 
% eight was flown and this is not enough to estimate parameters)
% see "importSimulationData"

function [meas_AccelerationKite0, meas_AccelerationKite1, meas_AccelerationKite2, ...
            meas_VelocityKite0, meas_VelocityKite1, meas_VelocityKite2, ...
            meas_PositionKite0, meas_PositionKite1, meas_PositionKite2, ...
            meas_MomentSumKite0, meas_MomentSumKite1, meas_MomentSumKite2, ...
            meas_AngularVelocityKite0, meas_AngularVelocityKite1, meas_AngularVelocityKite2, ...
            meas_EulerAnglesKite0,	meas_EulerAnglesKite1, meas_EulerAnglesKite2] ...
            = importData(file, ifPlotImportedData)

%'../../SimulationenBauer/20171108/Data.csv'
temp = csvread(file);

meas_AccelerationKite0 = temp(:,1);
meas_AccelerationKite1 = temp(:,2);
meas_AccelerationKite2 = temp(:,3);
meas_VelocityKite0 = temp(:,4);
meas_VelocityKite1 = temp(:,5);
meas_VelocityKite2 = temp(:,6);
meas_PositionKite0 = temp(:,7);
meas_PositionKite1 = temp(:,8);
meas_PositionKite2 = temp(:,9);
meas_MomentSumKite0 = temp(:,10);
meas_MomentSumKite1 = temp(:,11);
meas_MomentSumKite2 = temp(:,12);
meas_AngularVelocityKite0 = temp(:,13);
meas_AngularVelocityKite1 = temp(:,14);
meas_AngularVelocityKite2 = temp(:,15);
meas_EulerAnglesKite0 = temp(:,16);
meas_EulerAnglesKite1 = temp(:,17); 
meas_EulerAnglesKite2 = temp(:,18);

clear temp;

if ifPlotImportedData
    plotImportedData(meas_PositionKite0, meas_PositionKite1, meas_PositionKite2);
end

end