% Add blocks for impedance measurement in the created Simulink model.
% Run Usermain.m and create a simulink model first.
% Author(s): Ruiting Xu

%% Model Preparing
clear all; clc; close all; 
currentFile = mfilename('fullpath');
rootPattern = '+SimplusGT';
idx = strfind(currentFile, rootPattern)-length(rootPattern);
if ~isempty(idx)
    projectRoot = currentFile(1 : idx + length(rootPattern) - 1);
    targetPath = fullfile(projectRoot);
else
    error('Unable to locate project root directory automatically. Please specify the path manually.');
end

Folder_Name='fIMR_case';
load([targetPath,Folder_Name,'/','initial_data14.mat']);
Name_Model='IEEE14';
Name_LibFile = 'SimplusGT';
Tsample=Ts;

%% Add noise
Addnoise_flag=1;
Name_noisefile='Noise_raw';

%% Model Revising
% open the model
fprintf('==================================\n')
fprintf(['Opening the model: ',[Folder_Name,'/',Name_Model,'.slx']]); fprintf('\n')
fprintf('==================================\n')
open_system([Folder_Name,'/',Name_Model,'.slx']);

% some configs that might help improving efficiency
activeConfigObj = getActiveConfigSet(Name_Model);
set_param(activeConfigObj,'SaveTime','off');
set_param(activeConfigObj,'SaveOutput','off');
set_param(activeConfigObj,'SignalLogging','off');
set_param(activeConfigObj,'DSMLogging','off');
set_param(activeConfigObj,'ReturnWorkspaceOutputs','off');

% add purturbation and measurement blocks
Size_Theta = [40,30];
[Name_Thetablc, Name_Thetalabel] = SimplusGT.Simulink.SimAddTheta(Name_Model,Name_LibFile,Size_Theta,Pos_powergui);
Size_Perturb = [30,90];
ListPerturb=List_IBRbus;  %select measured buses
% ListPerturb=[1, 2, 3];  %select measured buses

% noise initialization
Noise_enable = 1;
Noise_raw = load(Name_noisefile);
Noise_raw = Noise_raw.Noise_raw;
nset = SimplusGT.Simulink.SimSetNoise(Noise_enable,Noise_raw,PowerFlow,ListPerturb,Fs);
Name_Perturb = SimplusGT.Simulink.SimAddPerturbation(Name_Model,Name_LibFile,Name_noisefile,Name_Thetalabel,...
    Size_Perturb,Pos_Bus,Pos_powergui,ListPerturb);
SimplusGT.Simulink.SimDisconnectApparatus_Bus(Name_Model,Name_Bus,Name_Apparatus,ListPerturb,ApparatusType);
SimplusGT.Simulink.SimConnectPerturbation(Name_Model,Name_Bus,Name_Apparatus,Name_Perturb,ListPerturb,ApparatusType);

% Limit the archive storage
Simulink.sdi.setArchiveRunLimit(15);

% also, save the entire workspace as "initial_data.mat" after running this script.
save([targetPath,Folder_Name,'/','initial_data14_1.mat']);