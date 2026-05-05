% Main function of frequency sweeping in simulation
% Authors: Ruiting Xu, Yue Zhu, etc.

% Before running, load initial parameters generated from usermain
clear all; clc; close all; 
load('initial_data14_1.mat'); 
 % selection of the model is set in CreateImpScanModel.m and FreqScan_main.m for flexiblity

% Step1: Create model for impedance sweeping
% SimplusGT.fIMR.CreateImpScanModel;

% Step2: Run Critical_fIMR.m for participation analysis
% SimplusGT.fIMR.FreqScan_main;

% Step 3: Plot scanned impedances and fIMR
load('Z_test_14.mat');
SimplusGT.fIMR.Z_test_plot;

% Step 4: Participation analysis;
if CaseStudy<=3
    List_Mode_sel = [30.7,39.15,64.81,76.43,124.41,631.408,1374.84];  
    No_modec = 4; 
    Appselc=10; %check List_IBRbus
else
    List_Mode_sel = [31.53,98.61];
    No_modec = 1; 
    Appselc=3; 
end
SimplusGT.fIMR.Critical_fIMR_scan;
