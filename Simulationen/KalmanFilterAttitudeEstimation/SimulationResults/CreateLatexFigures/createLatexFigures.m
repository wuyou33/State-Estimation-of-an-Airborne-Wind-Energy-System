% Creates using mat2tikz latex figures
% https://tomlankhorst.nl/matlab-to-latex-with-matlab2tikz/
% https://github.com/matlab2tikz/matlab2tikz
clc; clear; close all;

pathToMasterarbeitFolder = '../../../..';
pathToMat2Tikz = fullfile(pathToMasterarbeitFolder, 'Text/matlab2tikz-master/src');
pathFromLatexMain = '../Simulationen/KalmanFilterAttitudeEstimation/SimulationResults/CreateLatexFigures/';
pathToPresentationImagesDestinationFolder= fullfile(pathToMasterarbeitFolder, 'Präsentation/Abschlusspräsentation/Bilder');
addpath(pathToMat2Tikz)

exportForThesis = 0; % 0 is for the images for the presentation

%https://stackoverflow.com/questions/11212172/latex-fonts-in-matlab
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'DefaultAxesFontName', 'CMU Serif');

% paths to the *.fig files
files = {
        '../MagnetometerCalibration/magnetometerError.fig', ...
        %'../MagnetometerCalibration/MagnetometerHardIronCalibration.fig', ...
        %'../InputData/BauerDatenWithoutINCDEC.fig', ...
        %'../InputData/BauerDatenAttitude.fig', ...
        %'../../../KalmanFilterGleichstrommaschine/ResultPlots/100Hz.fig',...
        %'../../../KalmanFilterGleichstrommaschine/ResultPlots/ZoomIn.fig',...
        %'../../../KalmanFilterGleichstrommaschine/ResultPlots/ZoomOut.fig',...
        %'../ahrsEKF/angles.fig', ...
        %'../ahrsEKF/estBias.fig', ...
        %'../ahrsQuaternion/angles.fig', ...
        %'../ahrsMEKF/mekf.fig', ...
        %'../MahonyExplicitComplementaryFilter/angle.fig',...
        %'../INSMEKF/angles.fig', ...
        %'../INSMEKF/disableGPS4g.fig', ...
        %'../INSMEKF/gyroBias.fig', ...
        %'../INSMEKF/posVel.fig', ...
        %'../DifferentSampleTime/DST.fig', ...
        %'../AccelerationLimitedGPSModule/INSMEKF.fig',...
        %'../AccelerationLimitedGPSModule/Pos.fig',...
        %'../AccelerationLimitedGPSModule/GPSPos.fig',...
        %'../TetherSagCompensation/tetherSagCompensation.fig', ...
        %'../TetherSagCompensation/withoutTetherSagCompensationPosVel.fig', ...
        };

for i = 1 : length(files)
    file = files{i};
    fig = openfig(file);
    children = get(gcf, 'children');
    for i=1: length(children)
        child = children(i);
        if strcmp(child.Type, 'legend')
            %set(child, 'Visible', 'off');
        end
    end
    cleanfigure()
    [filepath,name,ext] = fileparts(file);
    texFileName = [name, '', '.tex'];
    texFile = fullfile(filepath, texFileName);
    % 'width', '\figW' not needed with standalone
    % use lualatex because pdf latex has problems with large data sets
    if exportForThesis
        relativePathFromMain = [pathFromLatexMain, filepath];
        matlab2tikz('figurehandle', fig, 'filename', texFile, ...
                'externalData', true, ...
                'tikzFileComment', ' !TeX program = lualatex',...
                'relativeDataPath', relativePathFromMain, 'width', '\figW');...
    else
        [filepath,name,ext] = fileparts(file);
        [~,folderName,~] = fileparts(filepath);
        name = fullfile(folderName, name);
        texFileName = [name, '', '.tex'];
        texFile = fullfile(pathToPresentationImagesDestinationFolder, texFileName);
        extraCode = ['\Large\newlength\figW\setlength\figW{12cm}'...
                    '\newlength\figH\setlength\figH{6cm}'];
        matlab2tikz('figurehandle', fig, 'filename', texFile, ...
                'externalData', true, ...
                'tikzFileComment', ' !TeX program = lualatex',...
                'width', '\figW', ...
                'height', '\figH', ...
                'standalone', true, 'extraCode', extraCode);
    end
    
    close(fig);
end
