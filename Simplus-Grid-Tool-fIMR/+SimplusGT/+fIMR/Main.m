% Author(s): Yitong Li, Yunjie Gu

%% 
% Notes:
%
% Please read "README.md" first before using the toolbox.
%
% Please use "Main_Customer.m" rather than this file for running the
% toolbox.

%%
fprintf('==================================\n')
fprintf('Start to run Simplus Grid Tool\n')
fprintf('==================================\n')

%% 
% ==================================================
% Change current path of matlab
% ==================================================
pathstr = mfilename('fullpath');        % Get the path of main.m
[pathstr,~,~]  = fileparts(pathstr);    % Get the path of Toolbox namespace
[pathstr,~,~]  = fileparts(pathstr);    % Get the of root namespace
[pathstr,~,~]  = fileparts(pathstr);    % Get the path of toolbox
cd(pathstr);                            % Change the current address

%%
% ==================================================
% Load customized data
% ==================================================
fprintf('Loading data, please wait a second...\n')

% ### Re-arrange basic settings
ListBasic = xlsread(UserData,'Basic');
Fs = ListBasic(1);
Ts = 1/Fs;              % (s), sampling period
Fbase = ListBasic(2);   % (Hz), base frequency
Sbase = ListBasic(3);   % (VA), base power
Vbase = ListBasic(4);   % (V), base voltage
Ibase = Sbase/Vbase;    % (A), base current
Zbase = Vbase/Ibase;    % (Ohm), base impedance
Ybase = 1/Zbase;        % (S), base admittance
Wbase = Fbase*2*pi;     % (rad/s), base angular frequency

% ### Re-arrange advanced settings
ListAdvance = xlsread(UserData,'Advance');
Flag_PowerFlowAlgorithm   	= ListAdvance(5);
Enable_CreateSimulinkModel	= ListAdvance(6);
Enable_PlotPole           	= ListAdvance(7);
Enable_PlotAdmittance     	= ListAdvance(8);
Enable_PrintOutput       	= ListAdvance(9);
Enable_Participation        = ListAdvance(10);

% ### Re-arrange the bus netlist
[ListBus,N_Bus] = SimplusGT.Toolbox.RearrangeListBus(UserData);

ListBus_=ListBus; % temparary added by Yue to save the original data

% ### Re-arrange the line netlist
[ListLine,N_Branch.N_Bus_] = SimplusGT.Toolbox.RearrangeListLine(UserData,ListBus);
DcAreaFlag = find(ListBus(:,12)==2);

% ### Re-arrange the apparatus netlist
[ApparatusBus,ApparatusType,Para,N_Apparatus] = SimplusGT.Toolbox.RearrangeListApparatus(UserData,Wbase,ListBus);
% The names of "ApparatusType" and "Para" can not be changed, because they
% will also be used in simulink model.

% Notes:
% No error checking if number of apparatuses is different from number of buses.

%%
% ==================================================
% Power flow analysis
% ==================================================

% ### Power flow analysis
fprintf('Doing the power flow analysis...\n')
if ~isempty(DcAreaFlag)
    Flag_PowerFlowAlgorithm = 1;
    fprintf(['Warning: Because the system has dc area(s), the Gauss-Seidel power flow method is always used.\n']);
end
switch Flag_PowerFlowAlgorithm
    case 1  % Gauss-Seidel 
        [PowerFlow,~,~,~,~,~,~,~] = SimplusGT.PowerFlow.PowerFlowGS(ListBus,ListLine,Wbase);
    case 2  % Newton-Raphson
       	[PowerFlow] = SimplusGT.PowerFlow.PowerFlowNR(ListBus,ListLine,Wbase);
    otherwise
        error(['Error: Wrong setting for power flow algorithm.']);
end
% Form of PowerFlow{i}: P, Q, V, xi, w
% P and Q are in load convention, i.e., the P and Q flowing from the bus to
% the apparatus.

% For printing later
ListPowerFlow = SimplusGT.PowerFlow.Rearrange(PowerFlow);

% Move load flow (PLi and QLi) to bus admittance matrix
[ListBus,ListLineNew,PowerFlowNew] = SimplusGT.PowerFlow.Load2SelfBranch(ListBus,ListLine,PowerFlow);

% For printting later
ListPowerFlowNew = SimplusGT.PowerFlow.Rearrange(PowerFlowNew);

%%
% ==================================================
% Descriptor state space model
% ==================================================

% ### Get the model of lines
fprintf('Getting the descriptor state space model of network lines...\n')

[YbusObj,YbusDSS,~] = SimplusGT.Toolbox.YbusCalcDss(ListBus,ListLineNew,Wbase);
[~,lsw] = size(YbusDSS.B);
ZbusObj = SimplusGT.ObjSwitchInOut(YbusObj,lsw);
[ZbusStateStr,ZbusInputStr,ZbusOutputStr] = ZbusObj.GetString(ZbusObj);

% ### Get the models of bus apparatuses
fprintf('Getting the descriptor state space model of bus apparatuses...\n')
for i = 1:N_Apparatus
    if length(ApparatusBus{i}) == 1
     	ApparatusPowerFlow{i} = PowerFlowNew{ApparatusBus{i}};
    elseif length(ApparatusBus{i}) == 2
        ApparatusPowerFlow{i} = [PowerFlowNew{ApparatusBus{i}(1)},PowerFlowNew{ApparatusBus{i}(2)}];
    else
        error(['Error']);
    end
    [GmObj_Cell{i},GmDSS_Cell{i},ApparatusPara{i},ApparatusEqui{i},ApparatusDiscreDamping{i},OtherInputs{i},ApparatusStateStr{i},ApparatusInputStr{i},ApparatusOutputStr{i}] = ...
        SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{i},ApparatusType{i},ApparatusPowerFlow{i},Para{i},Ts,ListBus);
    
    % The following data is not used in the script, but will be used in
    % simulations. Do not delete!
    x_e{i} = ApparatusEqui{i}{1};
    u_e{i} = ApparatusEqui{i}{2};
end

% ### Get the appended model of all apparatuses
fprintf('Getting the appended descriptor state space model of all apparatuses...\n')
GmObj = SimplusGT.Toolbox.ApparatusModelLink(GmObj_Cell);

% ### Get the model of whole system
fprintf('Getting the descriptor state space model of whole system...\n')
[GsysObj,GsysDSS,Port_v,Port_i,BusPort_v,BusPort_i] = ...
    SimplusGT.Toolbox.ConnectGmZbus(GmObj,ZbusObj,N_Bus);

% ### Whole-system admittance model
YsysObj = SimplusGT.ObjTruncate(GsysObj,Port_i,Port_v);
[~,YsysDSS] = YsysObj.GetDSS(YsysObj);

% SS string
GminSSObj = SimplusGT.ObjDss2Ss(GsysObj);
[GminStateStr,~,~]=GminSSObj.GetString(GminSSObj);

% ### Chech if the system is proper
fprintf('Checking if the whole system is proper:\n')
if isproper(GsysDSS)
    fprintf('Proper!\n');
    fprintf('Calculating the minimum realization of the system model for later use...\n')
    % GminSS = minreal(GsysDSS);
    GminSS = SimplusGT.dss2ss(GsysDSS);
    % This "minreal" function only changes the element sequence of state
    % vectors, but does not change the element sequence of input and output
    % vectors.
    InverseOn = 0;
else
    error('Error: System is improper, which has more zeros than poles.')
end
if SimplusGT.is_dss(GminSS)
    error(['Error: Minimum realization is in descriptor state space (dss) form.']);
end

%%
% ==================================================
% Create Simulink Model
% ==================================================
fprintf('\n')
fprintf('=================================\n')
fprintf('Simulink Model\n')
fprintf('=================================\n')

if N_Bus>=150
    Enable_CreateSimulinkModel = 0;
    fprintf('Warning: The system has more than 150 buses;\n')
    fprintf('         The simulink model can not be created because of the limited size of GUI.\n')
    fprintf('         The static and dynamic analysis will not be influenced.\n')
end

if Enable_CreateSimulinkModel == 1
    
    fprintf('Creating the simulink model automatically, please wait a second...\n')

    % Set the simulink model name
    Name_Model = 'mymodel_v1';

    % Close existing model with same name
    close_system(Name_Model,0);
    
    % Create the simulink model
    SimplusGT.Simulink.MainSimulink(Name_Model,ListBus,ListLineNew,ApparatusBus,ApparatusType,ListAdvance,PowerFlowNew);
    fprintf('Get the simulink model successfully! \n')
    fprintf('Please click the "run" button in the model to run it.\n')
    %fprintf('Warning: for later use of the simulink model, please "save as" a different name.\n')

else
    fprintf('Warning: The auto creation of simulink model is disabled.\n')
end

%%
% ==================================================
% Print results
% ==================================================

% ### Output the System
fprintf('\n')
fprintf('==================================\n')
fprintf('Print results\n')
fprintf('==================================\n')
fprintf('Whole system port model (system object form): GsysObj\n')
fprintf('Whole system port model (descriptor state space form): GsysDSS\n')
if Enable_PrintOutput
    [SysStateString,SysInputString,SysOutputString] = GsysObj.GetString(GsysObj);
    fprintf('Print ports of GsysDSS:\n')
    SimplusGT.Toolbox.PrintSysString(ApparatusBus,ApparatusType,GmObj_Cell,ZbusStateStr);
	fprintf('Print power flow result:\n')
    fprintf('The format below is "| bus | P | Q | V | angle | omega |". P and Q are in load convention.\n')
    ListPowerFlow
end

fprintf('Other models: \n')
fprintf('Whole system port model (state space form): GminSS\n')
fprintf('Whole system admittance model (system object form): YsysObj\n')
fprintf('Whole system admittance model (descriptor state space form): YsysDSS\n')

%%
% ==================================================
% Check Stability
% ==================================================

fprintf('\n')
fprintf('==================================\n')
fprintf('Check Stability\n')
fprintf('==================================\n')

fprintf('Calculatting pole/zero...\n')
pole_sys = pole(GsysDSS)/2/pi;
%fprintf('Checking if the system is stable:\n')
if isempty(find(real(pole_sys)>1e-8, 1))
    fprintf('Stable!\n');
else
    fprintf('Warning: Unstable!\n')
end

%%
% ==================================================
% Plot
% ==================================================

fprintf('\n')
fprintf('==================================\n')
fprintf('Plot Fundamentals\n')
fprintf('==================================\n')

% Initialize figure index
figure_n = 1000;

% Plot pole/zero map
if Enable_PlotPole
    fprintf('Plotting pole map...\n')
    figure_n = figure_n+1;
    figure(figure_n);
    subplot(1,2,1)
    scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Global pole map');
    
	subplot(1,2,2)
    scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    title('Zoomed pole map');
    axis([-80,20,-150,150]);
    plot([-80,0], [-80,0]*10, '--k','LineWidth',2,'Color','blue')
    plot([-80,0], [80,0]*10, '--k','LineWidth',2,'Color','blue')
    legend('mode','10% damping line')
    
    %SimplusGT.mtit('Pole Map');
else
    fprintf('Warning: The default plot of pole map is disabled.\n')
end

% omega_p = [logspace(-1,log10(25.5),5e2),logspace(log10(25.5),log10(27.5),50),logspace(log10(27.5),4,5e2)]*2*pi;
% omega_p = [logspace(-1,log10(90),5e2),logspace(log10(90),log10(92),1e2),logspace(log10(92),4,5e2)]*2*pi;
omega_p = [logspace(-1,4,1e3)]*2*pi;
omega_pn = [-flip(omega_p),omega_p];


Num_IBRbus=0; %number of the buses
for k = 1:N_Bus
    [k1,k2] = SimplusGT.CellFind(ApparatusBus,k);
    % if (0<=ApparatusType{k2} && ApparatusType{k2}<90) || ...
    %         (1000<=ApparatusType{k2} && ApparatusType{k2}<1090) || ...
    %         (2000<=ApparatusType{k2} && ApparatusType{k2}<2090)
    if ApparatusType{k}<30
        Num_IBRbus=Num_IBRbus+1;
        List_IBRbus(Num_IBRbus) = k;
    end
end

% Plot admittance
if Enable_PlotAdmittance
    fprintf('Plotting admittance spectrum, the calculation may take five minutes...')
    figure_n = figure_n+1;
    figure(figure_n);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:Num_IBRbus
        [k1,k2] = SimplusGT.CellFind(ApparatusBus,k);
        % Plot the active bus admittance only
        k_a=List_IBRbus(k);
        Yss{k_a}  = GminSS(BusPort_i{k_a},BusPort_v{k_a});
        Ysym{k_a} = SimplusGT.ss2sym(Yss{k_a});
        SimplusGT.bode_c(Ysym{k_a}(1,1),1j*omega_p,'PhaseOn',1,'LineWidth',0.5);
        CountLegend = CountLegend + 1;
        VecLegend{CountLegend} = ['POI-',num2str(k_a)];
    end
    legend(VecLegend);
    xlabel('Frequency (Hz)');
    % ylabel('Magnitude (pu)');
    SimplusGT.mtit('Admittance Spectrum');
else
    fprintf('Warning: The default plot of admittance spectrum is disabled.\n')
end

% Plot Ysys and Zapp
% OmegaPN = [-flip(omega_p),omega_p];
OmegaPN = omega_p;

for k=1:Num_IBRbus
    k_a=List_IBRbus(k);
    Ysysfun{k} = matlabFunction(Ysym{k_a});
    for i=1:length(OmegaPN)
        Yapp{k}(:,:,i) = freqresp(GmDSS_Cell{k_a}(1:2,1:2),1j*OmegaPN(i));
        Zapp_temp = inv(Yapp{k}(:,:,i));
        Zapp{k}(:,:,i) = Zapp_temp;
        Ysys_temp = Ysysfun{k}(1j*OmegaPN(i));
        Ysys{k}(:,:,i) = Ysys_temp;
        n_Ysys{k}(:,:,i) = norm(Ysys_temp);
        n_Zapp{k}(:,:,i) = norm(Zapp_temp);
        M{k}(:,:,i) = 1/(norm(Ysys_temp)*norm(Zapp_temp));
    end
end

axis_co=[1 300 1e-3 1e2];
figure_n=figure_n+1;
figure(figure_n)
for k=1:Num_IBRbus
    SimplusGT.plot_c(Zapp{k}(1,1,:),OmegaPN/2/pi,'PhaseOn',0,'PhaseShift',0,'LineWidth',0.5);
    % SimplusGT.bode_c(YcellSym{k}(1,1),1j*OmegaPN);
end
legend(VecLegend);
xlabel('Frequency (Hz)');
ylabel('Magnitude (pu)');
SimplusGT.mtit('Apparatus impedance');

figure_n=figure_n+1;
figure(figure_n)
for k=1:Num_IBRbus
    SimplusGT.plot_c(M{k}(1,1,:),OmegaPN/2/pi,'PhaseOn',0,'PhaseShift',0,'LineWidth',0.5);
    % SimplusGT.bode_c(YcellSym{k}(1,1),1j*OmegaPN);
end
legend(VecLegend);
xlabel('Frequency (Hz)');
ylabel('Magnitude (pu)');
SimplusGT.mtit('fIMR');
set(gcf,'unit','normalized','position',[0.2,0.2,0.2,0.16]);
set(gca,'position',[0.15,0.2,0.74,0.65]);
% set(gca,'XTickLabel',{' '});
set(gca,'fontsize',9,'fontname','Times New Roman');
set(gca,'YTick',[10^-4,10^-3,10^-2,10^-1,10^0,10^1,10^2]);
axis([1 1000 10^-2.5 10^0.9]);




