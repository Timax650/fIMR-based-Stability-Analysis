% The case studies for Frequency-spectra Impedance Margin Ratio based Stability Analysis of IBR-Penetrated Systems
% Authors: Ruiting Xu, Yue Zhu, etc.

% Notes: 
% * Add all folders and subfolders in the Matlab path before running.
% * My matlab version: R2024b
% * Certain Matlab add-ons may be required.

clear;
clc;
close all;

% choose the avaliable case studies: 1-2.
CaseStudy=1;
switch CaseStudy
    case 1; UserData = 'IEEE68.xlsx';  
    case 2; UserData = 'IEEE68_revised_para.xlsx';        
end
tic
SimplusGT.fIMR.Main();

% run CriticalIMRF.m;
toc 