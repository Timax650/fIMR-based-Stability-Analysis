% This function adds perturbation into the simulink model.

% Author(s): Ruiting Xu
function [Name_Perturb,Name_CurNoise,Name_VoltNoise] = SimAddPerturbation(Name_Model,Name_LibFile,Name_noisefile,...
                                                       Name_Thetalabel,Size_Perturb,Pos_Bus,Pos_Powergui,ListPerturb)

% Organize data
N_Perturb = length(ListPerturb);

Name_Perturb{1} = [];

% Mux, Scope, Toworkspace for Voltages
Name_MuxVolt = 'MuxVolt';
FullName_MuxVolt = [Name_Model '/' Name_MuxVolt];
add_block('simulink/Signal Routing/Mux',FullName_MuxVolt);
set_param(FullName_MuxVolt,'Inputs',num2str(N_Perturb));

Pos_MuxVolt = Pos_Powergui + [500+60,40*(-N_Perturb+1)-10];
size_MuxVolt=[5,40*N_Perturb];
set_param(FullName_MuxVolt,'position',[Pos_MuxVolt,Pos_MuxVolt+size_MuxVolt]);
set_param(gcb,'Orientation','right');

Name_VoltScope = 'VoltScope';
FullName_VoltScope = [Name_Model '/' Name_VoltScope];
add_block('simulink/Commonly Used Blocks/Scope',FullName_VoltScope);

Pos_VoltScope = Pos_Powergui + [500+100,40*(-N_Perturb+1)/2-5];
set_param(FullName_VoltScope,'position',[Pos_VoltScope,Pos_VoltScope+[30,30]]);
set_param(gcb,'Orientation','right');

Name_VoltTW = 'VoltTW';
FullName_VoltTW = [Name_Model '/' Name_VoltTW];
add_block('simulink/Sinks/To Workspace',FullName_VoltTW);
set_param(FullName_VoltTW,'VariableName','vdq_bus');
set_param(FullName_VoltTW,'SaveFormat','Array');

Pos_VoltTW = Pos_Powergui + [500+100,40*(-N_Perturb+1)/2+50];
set_param(FullName_VoltTW,'position',[Pos_VoltTW,Pos_VoltTW+[60,30]]);
set_param(gcb,'Orientation','right');

% Mux, Scope, Toworkspace for Apparatus Voltages
Name_MuxAVolt = 'MuxAVolt';
FullName_MuxAVolt = [Name_Model '/' Name_MuxAVolt];
add_block('simulink/Signal Routing/Mux',FullName_MuxAVolt);
set_param(FullName_MuxAVolt,'Inputs',num2str(N_Perturb));

Pos_MuxAVolt = Pos_Powergui + [700+60,40*(-N_Perturb+1)-10];
size_MuxAVolt=[5,40*N_Perturb];
set_param(FullName_MuxAVolt,'position',[Pos_MuxAVolt,Pos_MuxAVolt+size_MuxAVolt]);
set_param(gcb,'Orientation','right');

Name_AVoltScope = 'AVoltScope';
FullName_AVoltScope = [Name_Model '/' Name_AVoltScope];
add_block('simulink/Commonly Used Blocks/Scope',FullName_AVoltScope);

Pos_AVoltScope = Pos_Powergui + [700+100,40*(-N_Perturb+1)/2-5];
set_param(FullName_AVoltScope,'position',[Pos_AVoltScope,Pos_AVoltScope+[30,30]]);
set_param(gcb,'Orientation','right');

Name_AVoltTW = 'AVoltTW';
FullName_AVoltTW = [Name_Model '/' Name_AVoltTW];
add_block('simulink/Sinks/To Workspace',FullName_AVoltTW);
set_param(FullName_AVoltTW,'VariableName','Avdq_bus');
set_param(FullName_AVoltTW,'SaveFormat','Array');

Pos_AVoltTW = Pos_Powergui + [700+100,40*(-N_Perturb+1)/2+50];
set_param(FullName_AVoltTW,'position',[Pos_AVoltTW,Pos_AVoltTW+[60,30]]);
set_param(gcb,'Orientation','right');

% Mux, Scope, Toworkspace for Currents
Name_MuxCur = 'MuxCur';
FullName_MuxCur = [Name_Model '/' Name_MuxCur];
add_block('simulink/Signal Routing/Mux',FullName_MuxCur);
set_param(FullName_MuxCur,'Inputs',num2str(N_Perturb));

Pos_MuxVolt = Pos_Powergui + [900+60,40*(-N_Perturb+1)-10];
size_MuxVolt=[5,40*N_Perturb];
set_param(FullName_MuxCur,'position',[Pos_MuxVolt,Pos_MuxVolt+size_MuxVolt]);
set_param(gcb,'Orientation','right');

Name_CurScope = 'CurScope';
FullName_CurScope = [Name_Model '/' Name_CurScope];
add_block('simulink/Commonly Used Blocks/Scope',FullName_CurScope);

Pos_CurScope = Pos_Powergui + [900+100,40*(-N_Perturb+1)/2-5];
set_param(FullName_CurScope,'position',[Pos_CurScope,Pos_CurScope+[30,30]]);
set_param(gcb,'Orientation','right');

Name_CurTW = 'CurTW';
FullName_CurTW = [Name_Model '/' Name_CurTW];
add_block('simulink/Sinks/To Workspace',FullName_CurTW);
set_param(FullName_CurTW,'VariableName','idq_bus');
set_param(FullName_CurTW,'SaveFormat','Array');

Pos_CurTW = Pos_Powergui + [900+100,40*(-N_Perturb+1)/2+50];
set_param(FullName_CurTW,'position',[Pos_CurTW,Pos_CurTW+[60,30]]);
set_param(gcb,'Orientation','right');


% Perturbation blocks and labels
for i = 1:N_Perturb
    % Add Perturbation block
    Name_Perturb{i} = ['AC dq perturb' num2str(ListPerturb(i))];
    FullName_Perturb{i} = [Name_Model '/' Name_Perturb{i}];
    add_block([Name_LibFile '/AC dq perturb'],FullName_Perturb{i});

    Pos_Perturb{i} = Pos_Bus{ListPerturb(i)} + [-70,-20];
    set_param(FullName_Perturb{i},'position',[Pos_Perturb{i},Pos_Perturb{i}+Size_Perturb]);
    set_param(FullName_Perturb{i},'nset',['nset{',num2str(i),'}']);     
    set_param(gcb,'Orientation','right');

    % Add Perturbation label
    % theta
    Name_Label_thetafrom{i} = ['thetato' num2str(ListPerturb(i))];
    FullName_Label_thetafrom{i} = [Name_Model '/' Name_Label_thetafrom{i}];
    add_block('simulink/Signal Routing/From',FullName_Label_thetafrom{i});
    set_param(FullName_Label_thetafrom{i},'GotoTag',Name_Thetalabel);

    Pos_thetafrom{i} = Pos_Bus{ListPerturb(i)} + [-130,-40];
    set_param(FullName_Label_thetafrom{i},'position',[Pos_thetafrom{i},Pos_thetafrom{i}+[40,20]]);
    set_param(gcb,'Orientation','right');

    add_line(Name_Model,...
        {[Name_Label_thetafrom{i} '/1']}, {[Name_Perturb{i} '/1']},'autorouting','smart');

    % Perturbation voltage Goto
    Name_Label_Pvoltto{i} = ['Pvoltto' num2str(ListPerturb(i))];
    FullName_Label_Pvoltto{i} = [Name_Model '/' Name_Label_Pvoltto{i}];
    add_block('simulink/Signal Routing/Goto',FullName_Label_Pvoltto{i});
    set_param(FullName_Label_Pvoltto{i},'GotoTag',['Pvolt' num2str(ListPerturb(i))]);

    Pos_Pvoltto{i} = Pos_Bus{ListPerturb(i)} + [-15,-120];
    set_param(FullName_Label_Pvoltto{i},'position',[Pos_Pvoltto{i},Pos_Pvoltto{i}+[50,20]]);
    set_param(gcb,'Orientation','right');

    add_line(Name_Model,...
        {[Name_Perturb{i} '/1']}, {[Name_Label_Pvoltto{i} '/1']},'autorouting','smart');

    % Apparatus voltage Goto
    Name_Label_Avoltto{i} = ['Avoltto' num2str(ListPerturb(i))];
    FullName_Label_Avoltto{i} = [Name_Model '/' Name_Label_Avoltto{i}];
    add_block('simulink/Signal Routing/Goto',FullName_Label_Avoltto{i});
    set_param(FullName_Label_Avoltto{i},'GotoTag',['Avolt' num2str(ListPerturb(i))]);

    Pos_Avoltto{i} = Pos_Bus{ListPerturb(i)} + [-15,-80];
    set_param(FullName_Label_Avoltto{i},'position',[Pos_Avoltto{i},Pos_Avoltto{i}+[50,20]]);
    set_param(gcb,'Orientation','right');

    add_line(Name_Model,...
        {[Name_Perturb{i} '/2']}, {[Name_Label_Avoltto{i} '/1']},'autorouting','smart');

    % Perturbation current Goto
    Name_Label_Pcurto{i} = ['Pcurto' num2str(ListPerturb(i))];
    FullName_Label_Pcurto{i} = [Name_Model '/' Name_Label_Pcurto{i}];
    add_block('simulink/Signal Routing/Goto',FullName_Label_Pcurto{i});
    set_param(FullName_Label_Pcurto{i},'GotoTag',['Pcur' num2str(ListPerturb(i))]);

    Pos_Pcurto{i} = Pos_Bus{ListPerturb(i)} + [-15,-40];
    set_param(FullName_Label_Pcurto{i},'position',[Pos_Pcurto{i},Pos_Pcurto{i}+[45,20]]);
    set_param(gcb,'Orientation','right');

    add_line(Name_Model,...
        {[Name_Perturb{i} '/3']}, {[Name_Label_Pcurto{i} '/1']},'autorouting','smart');

    % gathering the labels with from
    %Perturbation voltage
    Name_Label_Pvoltfrom{i} = ['Pvoltfrom' num2str(ListPerturb(i))];
    FullName_Label_Pvoltfrom{i} = [Name_Model '/' Name_Label_Pvoltfrom{i}];
    add_block('simulink/Signal Routing/From',FullName_Label_Pvoltfrom{i});
    set_param(FullName_Label_Pvoltfrom{i},'GotoTag',['Pvolt' num2str(ListPerturb(i))]);

    Pos_Pvoltfrom{i} = Pos_Powergui + [500,40*(-N_Perturb+i)];
    set_param(FullName_Label_Pvoltfrom{i},'position',[Pos_Pvoltfrom{i},Pos_Pvoltfrom{i}+[50,20]]);
    set_param(gcb,'Orientation','right');

    %Apparatus voltage
    Name_Label_Avoltfrom{i} = ['Avoltfrom' num2str(ListPerturb(i))];
    FullName_Label_Avoltfrom{i} = [Name_Model '/' Name_Label_Avoltfrom{i}];
    add_block('simulink/Signal Routing/From',FullName_Label_Avoltfrom{i});
    set_param(FullName_Label_Avoltfrom{i},'GotoTag',['Avolt' num2str(ListPerturb(i))]);

    Pos_Avoltfrom{i} = Pos_Powergui + [700,40*(-N_Perturb+i)];
    set_param(FullName_Label_Avoltfrom{i},'position',[Pos_Avoltfrom{i},Pos_Avoltfrom{i}+[50,20]]);
    set_param(gcb,'Orientation','right');    

    %current
    Name_Label_Pcurfrom{i} = ['Pcurfrom' num2str(ListPerturb(i))];
    FullName_Label_Pcurfrom{i} = [Name_Model '/' Name_Label_Pcurfrom{i}];
    add_block('simulink/Signal Routing/From',FullName_Label_Pcurfrom{i});
    set_param(FullName_Label_Pcurfrom{i},'GotoTag',['Pcur' num2str(ListPerturb(i))]);

    Pos_Pcurfrom{i} = Pos_Powergui + [900,40*(-N_Perturb+i)];
    set_param(FullName_Label_Pcurfrom{i},'position',[Pos_Pcurfrom{i},Pos_Pcurfrom{i}+[45,20]]);
    set_param(gcb,'Orientation','right');

    % Connect Mux
    add_line(Name_Model,...
        {[Name_Label_Pvoltfrom{i} '/1']}, {[Name_MuxVolt '/' num2str(i)]},'autorouting','smart');
    add_line(Name_Model,...
        {[Name_Label_Avoltfrom{i} '/1']}, {[Name_MuxAVolt '/' num2str(i)]},'autorouting','smart');
    add_line(Name_Model,...
        {[Name_Label_Pcurfrom{i} '/1']}, {[Name_MuxCur '/' num2str(i)]},'autorouting','smart');
end

% Connect Scope and Toworkspace
add_line(Name_Model,...
    {[Name_MuxVolt '/1']}, {[Name_VoltScope '/1']},'autorouting','smart');
add_line(Name_Model,...
    {[Name_MuxAVolt '/1']}, {[Name_AVoltScope '/1']},'autorouting','smart');
add_line(Name_Model,...
    {[Name_MuxCur '/1']}, {[Name_CurScope '/1']},'autorouting','smart');
add_line(Name_Model,...
    {[Name_MuxVolt '/1']}, {[Name_VoltTW '/1']},'autorouting','smart');
add_line(Name_Model,...
    {[Name_MuxAVolt '/1']}, {[Name_AVoltTW '/1']},'autorouting','smart');
add_line(Name_Model,...
    {[Name_MuxCur '/1']}, {[Name_CurTW '/1']},'autorouting','smart');

end