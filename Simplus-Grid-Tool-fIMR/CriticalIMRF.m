N_Apparatus = evalin('base', 'N_Apparatus');
N_Bus = evalin('base', 'N_Bus');
ApparatusType = evalin('base', 'ApparatusType');
ApparatusBus = evalin('base', 'ApparatusBus');
ApparatusInputStr = evalin('base', 'ApparatusInputStr');
ApparatusOutputStr = evalin('base', 'ApparatusOutputStr');
ApparatusStateStr=evalin('base', 'ApparatusStateStr');
GminStateStr=evalin('base', 'GminStateStr');
GminSS = evalin('base', 'GminSS');
GmDSS_Cell = evalin('base', 'GmDSS_Cell');
GsysDSS = evalin('base', 'GsysDSS');
Para = evalin('base', 'Para');
ApparatusPowerFlow = evalin('base', 'ApparatusPowerFlow');
Ts = evalin('base', 'Ts');
ListBus = evalin('base', 'ListBus');

GmObj=evalin('base', 'GmObj');
YbusObj=evalin('base', 'YbusObj');
Port_i=evalin('base', 'Port_i');
Port_v=evalin('base', 'Port_v');


%% enter modes and apparatus of discussion here
 % discussed modes for generating the heat map
kc = searchmode([26.14,30.37,90.77,143.9],OmegaPN/2/pi);    %68-bus, original
% kc = searchmode([26.14,43.9,109.154,143.9],OmegaPN/2/pi);    %68-bus, retuned
 % the discussed mode and apparatus for sensitivity analysis
No_mode = 3; %the second mode(30.37Hz)
Appsel=2; %2,4,5 → 26,28,29

k_wc = kc;
wc = OmegaPN(k_wc)/2/pi;

%% config C-IMR: Floating bus with maximum IMR value, SG with no effect on heat map, and IBR has the critical effect.
j=1; Num_app=0; 
clear CfIMR;
for k=1:N_Bus
    if ApparatusType{k}<30 || ApparatusType{k}>=40
        CfIMR(j).device = k;
        CfIMR(j).value = log10(100);
        CfIMR(j).mode = 0;
        j=j+1;
    end
    if ApparatusType{k}~=100
        Num_app=Num_app+1;
        app_list(Num_app)=k;
    end
end

for k=1:Num_IBRbus
    M_sq(k,:)=squeeze(M{k});
end

%% sweep the mode
for modei=1:length(wc)
    j_a=0;
    for j=1:length(CfIMR)
        k=CfIMR(j).device;
        if any(k == List_IBRbus)
            j_a=j_a+1;
            fIMR = M_sq(j_a,k_wc(modei));
            fIMR_o=fIMR;
            if fIMR<0.01
                fIMR = log10(0.01);
            else
                fIMR = log10(fIMR);
            end

            if fIMR<CfIMR(j).value
                CfIMR(j).value=fIMR;   
                CfIMR(j).mode = wc(modei);
                CfIMR(j).value_orig=fIMR_o;
            end
        end
    end
end


%% Output Ysys and fIMR at the discussed modes
Ysys_md_val=zeros(Num_IBRbus,length(wc));
Ysys_md_abs=zeros(Num_IBRbus,length(wc));
IMR_md_val=zeros(Num_IBRbus,length(wc));
contr=zeros(Num_IBRbus,length(wc));
realcontr_pu=zeros(Num_IBRbus,length(wc));
imagcontr_pu=zeros(Num_IBRbus,length(wc));
contrRealSum=zeros(1,length(wc));
contrImagSum=zeros(1,length(wc));
for i=1:Num_IBRbus
    for j=1:length(wc)
        Ysys_md_val(i,j)=Ysys{i}(1,1,k_wc(j));
        Ysys_md_abs(i,j)=abs(Ysys_md_val(i,j));
        IMR_md_val(i,j)=M_sq(i,k_wc(j));
        contr(i,j) = Ysys{i}(1,1,k_wc(j))*Zapp{i}(1,1,k_wc(j)) + Ysys{i}(1,2,k_wc(j))*Zapp{i}(1,2,k_wc(j)) +...
                     Ysys{i}(2,1,k_wc(j))*Zapp{i}(2,1,k_wc(j)) + Ysys{i}(2,2,k_wc(j))*Zapp{i}(2,2,k_wc(j));  
        contrRealSum(j) = contrRealSum(j) + abs(real(contr(i,j)));
        contrImagSum(j) = contrImagSum(j) + abs(imag(contr(i,j)));
    end
end
 %normalize
for j=1:length(wc)
    realcontr_pu(:,j)=real(contr(:,j))/contrRealSum(j);
    imagcontr_pu(:,j)=imag(contr(:,j))/contrImagSum(j);
end

%% plot
 %heatmap
Main_K_Plot(CfIMR,1);
title('Small-Signal System Strength Heatmap');

 %contribution factor
No_pole_plot=No_mode;
barwidth=0.7;
figure_n2=figure_n+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.1]);
b1=bar(realcontr_pu(:,No_pole_plot),barwidth);
b1.EdgeColor = 'none';
axis([0 15 -1 1]);
set(gca,'XTicklabel',List_IBRbus_plot);
set(gca,'position',[0.15,0.25,0.8,0.6]);
set(gca,'fontsize',9,'fontname','Times New Roman');
title('Contribution factor (Real)');

figure_n2=figure_n2+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.1]);
b2=bar(imagcontr_pu(:,No_pole_plot),barwidth);
b2.EdgeColor = 'none';
axis([0 15 -1 1]);
for i=1:length(List_IBRbus)
    List_IBRbus_plot{i}=num2str(List_IBRbus(i));
end
set(gca,'XTicklabel',List_IBRbus_plot);
set(gca,'position',[0.15,0.25,0.8,0.6]);
% set(gca,'fontsize',Fsize,'fontname','Times New Roman');
set(gca,'fontsize',9,'fontname','Times New Roman');
title('Contribution factor (Imaginary)');

%% sensitivity analysis
pole_sel = 1j*wc(No_mode);  %in Hz
AppSelect = List_IBRbus; 
for k=1:Num_IBRbus
    Ysys_sel{k} = Ysys{k}(:,:,k_wc(No_mode));
end
MdLayer3 = SimplusGT.fIMR.MdLayer3f(Ysys_sel,pole_sel,ApparatusType,...
        AppSelect,Para,ApparatusPowerFlow,Ts,ApparatusBus,ListBus);  

figure_n2=figure_n2+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.2,0.17]);
hold on
for parasel=1:length(MdLayer3(Appsel).Result)
    Imp_Para_sens(parasel)=MdLayer3(Appsel).Result(parasel).DeltaZ;
    Det_Para_sens(parasel)=MdLayer3(Appsel).Result(parasel).Ddet;
    plot([0,real(Det_Para_sens(parasel))],[0,imag(Det_Para_sens(parasel))]);
end
hold off
set(gca,'position',[0.15,0.2,0.7,0.7]);
axis([-12 12 -7 5])
set(gca,'XTick',-10:2:10);
set(gca,'YTick',-10:2:10);
set(gca,'fontsize',9,'fontname','Times New Roman');
title('Determinant-parameter sensitivity');
xlabel('real')
ylabel('imaginary');

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