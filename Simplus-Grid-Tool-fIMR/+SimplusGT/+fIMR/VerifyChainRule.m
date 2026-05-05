% Verification of the chain rule
% Run UserMain and Critical_fIMR first
%% Create state space mode
% cd(pathstr);
fprintf('Recreate Objs...\n')
clear GmObj_Cell YbusObj ZbusObj
for i = 1:N_Apparatus
    [GmObj_Cell{i},~,~,~,~,~,~,~,~] = ...
        SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{i},ApparatusType{i},ApparatusPowerFlow{i},Para{i},Ts,ListBus);
end
[YbusObj,~,~] = SimplusGT.Toolbox.YbusCalcDss(ListBus,ListLineNew,Wbase);
[~,lsw] = size(YbusDSS.B);
ZbusObj = SimplusGT.ObjSwitchInOut(YbusObj,lsw);


%% Verifycation
wvc = fc(No_modev)*2*pi;
Appsel_v = find(app_list == List_IBRbus(Appselv));   %Select in app_list

%calculate Zsum and Zdet
YbusSS = YbusDSS;
ListBus_= ListBus(:,1);
clear ZappDSS_Cell
n=0;
for k_a = app_list
    n=n+1;
    ZappDSS_Cell{n} = inv(GmDSS_Cell{k_a}(1:2,1:2));
end
ZappSS = append(ZappDSS_Cell{:});
Zapp_f = freqresp(ZappSS,wvc);
Ybus_f = freqresp(YbusSS,wvc);
Ybus_kron_f = kronZN(Ybus_f,app_list,ListBus_);
Zbus_kron_f = inv(Ybus_kron_f);
Zsum_f = Zapp_f+Zbus_kron_f; 
Ysys_f = inv(Zsum_f);
Zdet_f = det(Zsum_f);

mode_all_orig = pole_sys;
mode_pick_orig = pickmode(mode_all_orig,mode_pickv);


%Revise parameter
ParamName = fieldnames(Para{app_list(Appsel_v)});
ParaSel = getfield(Para{app_list(Appsel_v)},ParamName{Paraselv}); % extract the parameter
delta_para = 1e-2*(abs(ParaSel));
ParaPerturb = ParaSel + delta_para ; % add perturabation
ParaNew = setfield(Para{app_list(Appsel_v)}, ParamName{Paraselv}, ParaPerturb); % update the parameter
[GmObj_Cellk_New,GmDSS_Cellk_New,~,~,~,~,~,~,~] ...
    = SimplusGT.Toolbox.ApparatusModelCreate(app_list(Appsel_v),ApparatusType{app_list(Appsel_v)},...
    ApparatusPowerFlow{app_list(Appsel_v)},ParaNew,Ts,ListBus);

%calculate new Zapp and Zsum
ZappDSS_Cell_New=ZappDSS_Cell;
ZappDSS_Cell_New{Appsel_v}=inv(GmDSS_Cellk_New(1:2,1:2));
ZappSS_New = append(ZappDSS_Cell_New{:});
Zapp_f_New = freqresp(ZappSS_New,wvc);
Zsum_f_New = Zapp_f_New+Zbus_kron_f;

%Predicted result with sensitivity
delta_Zapp = Zapp_f_New(2*Appsel_v-1:2*Appsel_v,2*Appsel_v-1:2*Appsel_v)-Zapp_f(2*Appsel_v-1:2*Appsel_v,2*Appsel_v-1:2*Appsel_v);
Ysys_f_ = Ysys_f(2*Appsel_v-1:2*Appsel_v,2*Appsel_v-1:2*Appsel_v);
delta_Zdet_pred = Ysys_f_(1,1)*delta_Zapp(1,1) + Ysys_f_(1,2)*delta_Zapp(2,1) +...
                  Ysys_f_(2,1)*delta_Zapp(1,2) + Ysys_f_(2,2)*delta_Zapp(2,2);
Zdet_f_pred = delta_Zdet_pred*Zdet_f + Zdet_f;
mode_pred = mode_pick_orig + delta_Zdet_pred*(mode_pick_orig*2*pi-1i*wvc)/2/pi;

%calculate new Zdet and poles
Zdet_f_New = det(Zsum_f_New);

GmObj_Cell_New=GmObj_Cell;
GmObj_Cell_New{app_list(Appsel_v)}=GmObj_Cellk_New;
GmObj_New = SimplusGT.Toolbox.ApparatusModelLink(GmObj_Cell_New);
[~,GsysDSS_New,~,~,~,~] = ...
  SimplusGT.Toolbox.ConnectGmZbus(GmObj_New,ZbusObj,N_Bus);

mode_all_New = pole(GsysDSS_New)/2/pi;
mode_pick_New = pickmode(mode_all_New,mode_pickv);

%check establishment of determinant-mode sensitivity
G_orig=Zdet_f/(1i*wvc-mode_pick_orig*2*pi)/2/pi;  G_New=Zdet_f_New/(1i*wvc-mode_pick_New*2*pi)/2/pi;
delt_G=abs((G_New-G_orig)/G_orig);
delt_mode=abs(2*pi*(mode_pick_New-mode_pick_orig)/(1i*wvc-mode_pick_orig*2*pi));
% checkZdet=(Zdet_f_New-Zdet_f)-((1i*wvc-mode_pick_orig*2*pi)*(G_New-G_orig)-G_orig*2*pi*(mode_pick_New-mode_pick_orig));
fprintf(['ΔG/G = ' num2str(delt_G) '\n'])
fprintf(['Δλ/(s-λ) = ' num2str(delt_mode) '\n'])
if delt_G<0.1*delt_mode
    fprintf(['Determinant-mode sensitivity is accurate\n'])
else
    fprintf(['Determinant-mode sensitivity is not accurate enough\n'])
end


%plot result
 %Zdet
figure_n3=figure_n2+1;
figure(figure_n3); clf;
quiverB=quiver(real(Zdet_f),imag(Zdet_f),real(Zdet_f_New-Zdet_f),imag(Zdet_f_New-Zdet_f),'g-','AutoScale','off');
hold on
quiverC=quiver(real(Zdet_f),imag(Zdet_f),real(Zdet_f_pred-Zdet_f),imag(Zdet_f_pred-Zdet_f),'r--','AutoScale','off');
pointA=plot(Zdet_f,'b*');
pointB=plot(Zdet_f_New,'g*');
pointC=plot(Zdet_f_pred,'ro');
hold off
set(gcf,'unit','normalized','position',[0.2,0.2,0.15,0.275]);
set(gca,'position',[0.18,0.2,0.74,0.7]);
set(gca,'fontsize',14,'fontname','Times New Roman');
% x_start=-5.6e-68;
% y_start=-4.4e-68;
% sizeax=0.99e-68;
% axis([x_start x_start+sizeax y_start y_start+sizeax]);

% plot settings
Blue=[84, 139, 212]/255; Green=[84, 212, 139]/255; Red=[255, 0, 0]/255;
Linewidth_all=1.5; Markersize_all=10; MaxHeadSize_all=0.5;
set(pointA, 'Color', Blue);
set(pointA, 'LineWidth', Linewidth_all);
set(pointA, 'Markersize', Markersize_all);  
set(pointB, 'Color', Green); 
set(pointB, 'LineWidth', Linewidth_all); 
set(pointB, 'Markersize', Markersize_all);  
set(pointC, 'Color', Red); 
set(pointC, 'LineWidth', Linewidth_all); 
set(pointC, 'Markersize', Markersize_all);  
set(quiverB, 'Color', Green);
set(quiverB, 'LineWidth', Linewidth_all);
set(quiverB, 'MaxHeadSize', MaxHeadSize_all); 
set(quiverC, 'Color', Red); 
set(quiverC, 'LineWidth', 1.5); 
set(quiverC, 'MaxHeadSize', MaxHeadSize_all); 

 %mode
figure_n3=figure_n3+1;
figure(figure_n3); clf;
quiverB=quiver(real(mode_pick_orig),imag(mode_pick_orig),real(mode_pick_New-mode_pick_orig),imag(mode_pick_New-mode_pick_orig),'g-','AutoScale','off');
hold on
quiverC=quiver(real(mode_pick_orig),imag(mode_pick_orig),real(mode_pred-mode_pick_orig),imag(mode_pred-mode_pick_orig),'r--','AutoScale','off');
pointA=plot(mode_pick_orig,'b*');
pointB=plot(mode_pick_New,'g*');
pointC=plot(mode_pred,'ro');

hold off
set(gcf,'unit','normalized','position',[0.2,0.2,0.15,0.275]);
set(gca,'position',[0.18,0.2,0.74,0.7]);
set(gca,'fontsize',14,'fontname','Times New Roman');
% x_start=-1.86;
% y_start=30.79;
% sizeax=0.25;
% axis([x_start x_start+sizeax y_start y_start+sizeax]);

% plot settings
Blue=[84, 139, 212]/255; Green=[84, 212, 139]/255; Red=[255, 0, 0]/255;
Linewidth_all=1.5; Markersize_all=10; MaxHeadSize_all=0.5;
set(pointA, 'Color', Blue);
set(pointA, 'LineWidth', Linewidth_all);
set(pointA, 'Markersize', Markersize_all);  
set(pointB, 'Color', Green); 
set(pointB, 'LineWidth', Linewidth_all); 
set(pointB, 'Markersize', Markersize_all);  
set(pointC, 'Color', Red); 
set(pointC, 'LineWidth', Linewidth_all); 
set(pointC, 'Markersize', Markersize_all);  
set(quiverB, 'Color', Green);
set(quiverB, 'LineWidth', Linewidth_all);
set(quiverB, 'MaxHeadSize', MaxHeadSize_all); 
set(quiverC, 'Color', Red); 
set(quiverC, 'LineWidth', 1.5); 
set(quiverC, 'MaxHeadSize', MaxHeadSize_all); 

%  retrun to the original parameters
clear GmObj_Cell
for i = 1:N_Apparatus
    [GmObj_Cell{i},~,~,~,~,~,~,~,~] = ...
        SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{i},ApparatusType{i},ApparatusPowerFlow{i},Para{i},Ts,ListBus);
end



function kc=searchmode(wc,frequency)
kc=zeros(1,length(wc));
for i=1:length(wc)
    for j=1:length(frequency)
        if j>1
            if frequency(j-1)<wc(i) && frequency(j)>=wc(i)
                kc(i)=j;
            end
        end
    end
end
end

function ZN_New = kronZN(ZN,app_list,ListBus)
for i=1:length(app_list)
    keep_nodes(2*i-1)=2*app_list(i)-1;
    keep_nodes(2*i)=2*app_list(i);
end
for i=1:length(ListBus)
    all_nodes(2*i-1)=2*ListBus(i)-1;
    all_nodes(2*i)=2*ListBus(i);
end

elim_nodes = setdiff(all_nodes, keep_nodes);

Zrr = ZN(keep_nodes, keep_nodes);
Zre = ZN(keep_nodes, elim_nodes);
Zer = ZN(elim_nodes, keep_nodes);
Zee = ZN(elim_nodes, elim_nodes);

ZN_New = Zrr - Zre * (Zee \ Zer);
end


function mode_pick = pickmode(mode_all,mode_Hz)
    distances = abs(mode_all - mode_Hz);
    [~, index] = min(distances(:));
    mode_pick = mode_all(index);
end