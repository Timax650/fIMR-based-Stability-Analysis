% The case studies for Frequency-spectra Impedance Margin Ratio based Stability Analysis of IBR-Penetrated Systems
% Authors: Ruiting Xu, Yue Zhu, etc.

% Notes: 
% * Add all folders and subfolders in the Matlab path before running.
% * My matlab version: R2024b
% * Certain Matlab add-ons may be required.

clear all;
clc;
close all;

%% choose the avaliable case studies: 1-2.
CaseStudy=4;
switch CaseStudy
    case 1; UserData = 'IEEE68_GFM_original.xlsx';  
    case 2; UserData = 'IEEE68_GFM_detuned.xlsx';  
    case 3; UserData = 'IEEE68_GFM_improved.xlsx';     
    case 4; UserData = 'IEEE14_GFM_detuned.xlsx';    
    case 5; UserData = 'IEEE14_GFM4_m1.xlsx';
    case 6; UserData = 'IEEE14_GFM4_m.xlsx';
    case 7; UserData = 'IEEE68_GFM_m.xlsx';
    case 8; UserData = 'IEEE68_GFM_m1.xlsx';
end

%% Step1: Run Main.m to calculate Ysys and fIMR
tic
SimplusGT.fIMR.Main();
toc
newcase=1;
if newcase
    Wbase_sim=Wbase;
    Folder_Name='fIMR_case';
    % save([Folder_Name,'/','initial_data68.mat']);
    save([Folder_Name,'/','initial_data14.mat']);
end
  % Data calculated from main.m is the foundation of all analysis and models, load and save these data for convenience 

%% Step2: Run Critical_fIMR.m for participation analysis
load('initial_data14.mat');
if CaseStudy<=3
    List_Mode_sel = [30.7,39.15,64.81,76.43,124.41,631.408,1374.84];  
    No_modec = 4; 
    Appselc=10; %check List_IBRbus
else
    List_Mode_sel = [31.53,98.61];
    No_modec = 1; 
    Appselc=3; 
end
SimplusGT.fIMR.Critical_fIMR();
 %Run Critical_IMR.m for comparative verification
% if CaseStudy<=3
%     Modesel_=-1.666+30.844i;
% else
%     Modesel_=-0.584+31.73i;
% end
% SimplusGT.fIMR.Critical_IMR();
 % Save the workspace
save([Folder_Name,'/','initial_data14_2.mat']);

%% Step 3: Run VerifyChainRule.m for verification, the manuscript illustrates results of 68-bus system
load('initial_data14_2.mat');
No_modev=No_modec; Appselv=Appselc;
Paraselv = 8;   %Select any parameter which can change Zapp
if CaseStudy<=3
    mode_pickv = -1.667+30.84i; %in Hz, select the mode you want to verify, which should be aligned with the one you select in Critical_fIMR
    % mode_pick = -0.90+76.45i;
else
    mode_pickv = -0.584+31.73i;
end

SimplusGT.fIMR.VerifyChainRule();
