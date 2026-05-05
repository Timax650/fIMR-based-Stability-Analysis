% Script for frequency sweeping.
% Open a model created with CreateImpScanModel.m.

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

Folder_Name='fIMR_case'; Model_Name='IEEE14_1';
load([targetPath,Folder_Name,'/','initial_data14_1.mat']);  %save the workspace the first time you run Usermain.m and create a simulink model
t_init=0.1; %Time for simulating the steaty state
Tsample=Ts; 

%set frequency of sweeping
% frequencys=freq_gen(75,200,0.05);
fstart=400; fend=500; fdevmax=0.1; fdevmin=fdevmax/100; 
% regulate magnitude of the reponse current around 0.05pu
i_idea=0.05;

clear imp_meas

fprintf('==================================\n')
fprintf(['Opening the model: ',[Folder_Name,'/',Model_Name,'.slx']]); fprintf('\n')
fprintf('==================================\n')
open_system([Folder_Name,'/',Model_Name,'.slx']);
activeConfigObj = getActiveConfigSet(Model_Name);

N_Perturb = length(ListPerturb);
for N=1:N_Perturb
    dis_module_ac{N}=['/AC dq perturb' num2str(ListPerturb(N))];
end

fprintf('==================================\n')
fprintf('Run the model for the initial point\n')
fprintf('==================================\n')
for N=1:N_Perturb
    perturbation_set(Model_Name,dis_module_ac{N},0,1,'d');
    perturbation_set(Model_Name,dis_module_ac{N},0,1,'q');
end
set_param(activeConfigObj,'StopTime',num2str(t_init));
set_param(activeConfigObj,'SaveFinalState','on');
set_param(activeConfigObj,'SaveOperatingPoint','off');
set_param(activeConfigObj,'LoadInitialState','off');
set_param(activeConfigObj,'FinalStateName','x_init');
sim(Model_Name);

fprintf('==================================\n')
fprintf('Frequency Scanning...\n')
fprintf('==================================\n')
set_param(activeConfigObj,'SaveFinalState','off');
set_param(activeConfigObj,'LoadInitialState','on');
set_param(activeConfigObj,'InitialState','x_init');

Vcplx1=zeros(2,N_Perturb); Icplx1=zeros(2,N_Perturb);
Vcplx2=zeros(2,N_Perturb); Icplx2=zeros(2,N_Perturb);
Y_last=5*ones(2,N_Perturb); 

n=1; frequencys(1)=fstart; fdev=fdevmin; var_YZ=1e-9;  %initialization sweeping steps
while frequencys(n)<=fend
    fprintf(['Time: ' num2str(n) '      Frequency: ' num2str(frequencys(n)) '\n'])
    for N=1:N_Perturb
        fprintf(['Bus No.: ' num2str(N) '\n'])
        %
        [t_scan,ts]=ts_gen(frequencys(n),Ts); 
        set_param(activeConfigObj,'StopTime',num2str(t_scan));
        
        %d axis injection
        dq_axis = 'd';
        amplitude=FreqScan.Amplitude_choose(frequencys(n),Y_last(1,N),i_idea); amplitudes_ob(1,N,n)=amplitude; 
        perturbation_set(Model_Name,dis_module_ac{N},amplitude,frequencys(n),dq_axis);
        for N2=1:N_Perturb
            if N2~=N
                perturbation_set(Model_Name,dis_module_ac{N2},0,frequencys(n),dq_axis);
            end
        end
        sim(Model_Name);
        %FFT       
        Vdq_t1=vdq_bus(ceil((t_scan-ts)/Tsample)+1:end,2*N-1:2*N);
        AVdq_t1=Avdq_bus(ceil((t_scan-ts)/Tsample)+1:end,2*N-1:2*N);
        Idq_t1=idq_bus(ceil((t_scan-ts)/Tsample)+1:end,2*N-1:2*N);
        Vdq1(:,N)=FreqScan.FFT_analyze(Vdq_t1,frequencys(n),ts);
        AVdq1(:,N)=FreqScan.FFT_analyze(AVdq_t1,frequencys(n),ts);
        Idq1(:,N)=FreqScan.FFT_analyze(Idq_t1,frequencys(n),ts);

        %q axis injection
        dq_axis = 'q';
        amplitude=FreqScan.Amplitude_choose(frequencys(n),Y_last(2,N),i_idea); amplitudes_ob(2,N,n)=amplitude;
        perturbation_set(Model_Name,dis_module_ac{N},amplitude,frequencys(n),dq_axis);
        for N2=1:N_Perturb
            if N2~=N
                perturbation_set(Model_Name,dis_module_ac{N2},0,frequencys(n),dq_axis);
            end
        end        
        sim(Model_Name);
        %FFT
        Vdq_t2=vdq_bus(ceil((t_scan-ts)/Tsample)+1:end,2*N-1:2*N);
        AVdq_t2=Avdq_bus(ceil((t_scan-ts)/Tsample)+1:end,2*N-1:2*N);
        Idq_t2=idq_bus(ceil((t_scan-ts)/Tsample)+1:end,2*N-1:2*N);
        Vdq2(:,N)=FreqScan.FFT_analyze(Vdq_t2,frequencys(n),ts);
        AVdq2(:,N)=FreqScan.FFT_analyze(AVdq_t2,frequencys(n),ts);
        Idq2(:,N)=FreqScan.FFT_analyze(Idq_t2,frequencys(n),ts);

        %calculate impedance and save
        Y(:,:,N)=[Idq1(1,N),Idq2(1,N);Idq1(2,N),Idq2(2,N)]*...
            ([Vdq1(1,N),Vdq2(1,N);Vdq1(2,N),Vdq2(2,N)]\eye(2));
        Y_last(1,N)=Y(1,1,N);   Y_last(2,N)=Y(2,2,N);
        ZA(:,:,N)=[AVdq1(1,N),AVdq2(1,N);AVdq1(2,N),AVdq2(2,N)]*...
            ([Idq1(1,N),Idq2(1,N);Idq1(2,N),Idq2(2,N)]\eye(2));
        imp_meas(n).frequency=frequencys(n);
        imp_meas(n).Ysys=Y;
        imp_meas(n).ZA=ZA;

        save([targetPath,Folder_Name,'/','Z_test.mat'],'imp_meas');
    end
    if n > 1
        var_YZ=var_YZ_Calc(imp_meas(n).Ysys,imp_meas(n-1).Ysys,frequencys(n));
        fdev=max(min(freq_step_gen(fdev,var_YZ),fdevmax),fdevmin);
    end
    fdev_record(n)=fdev;
    n=n+1;
    frequencys(n)=10^( log10(frequencys(n-1))+fdev );
end

function perturbation_set(Model_Name,dis_module_ac,amplitude,frequency,dq_axis)
frequency_dis=frequency;
switch dq_axis
    case 'd'
        set_param([Model_Name,dis_module_ac],'fr',num2str(frequency_dis));
        set_param([Model_Name,dis_module_ac],'amp',['[',num2str(amplitude),',0]']);
    case 'q'
        set_param([Model_Name,dis_module_ac],'fr',num2str(frequency_dis));
        set_param([Model_Name,dis_module_ac],'amp',['[0,',num2str(amplitude),']']);
end
end

function f=freq_gen(fstart,fend,fdev)
log_f=log10(fstart):fdev:log10(fend);
f=zeros(1,length(log_f));

for k=1:length(log_f)
    if log_f(k)>1
        f(k)=floor(10^(log_f(k)));
    else
        f(k)=floor(10*10^(log_f(k)))/10;
    end

    % f(k)=ceil(10^(log_f(k)));
end

i=2;
while i<=length(f)
    if f(i)==f(i-1)
        f(i)=[];
        i=i-1;
    end
    i=i+1;
end
end

function [t_scan,ts]=ts_gen(frequency,Ts)
% if frequency<=2
%     t_scan=2/frequency; %É¨ÆµÊ±³¤
% elseif frequency<=4
%     t_scan=4/frequency;
% elseif frequency<=10
%     t_scan=6/frequency;
% else
%     t_scan=60/frequency;
% end
if frequency<=10
    t_scan=10/frequency;
else
    t_scan=1/10;
end
ts=2/frequency;
end

function var_YZ=var_YZ_Calc(Y,Ylast,freq_dis)
if freq_dis<200
    var_YZ=max(abs((Y-Ylast)./Ylast),[],"all");
else
    diagonals = zeros(2, size(Y,3));
    for k = 1:size(Y,3)
        diagonals(:, k) = diag(abs((Y(:,:,k)-Ylast(:,:,k))./Ylast(:,:,k)));
    end
    var_YZ=max(diagonals,[],"all"); %In high frequency band, coupling terms are very small
end
end


function freq_step=freq_step_gen(freq_step_last,var_YZ)
freq_step = freq_step_last * (0.4/var_YZ); %controlled impedance variation between each point: 40%
end