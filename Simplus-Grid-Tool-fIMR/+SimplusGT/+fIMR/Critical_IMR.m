% participation and sensitivity analysis based on the IMR method, for comparison 
% see papers:
% Y. Zhu, T. C. Green, et al, “Impedance Margin Ratio: a New Metric for Small-Signal System Strength,” IEEE Transactions on Power Systems, 2024, doi: 10.1109/TPWRS.2024.3371231.
%

%select mode and apparatus for participation analyais
A=GminSS.A;
[~,D]=eig(A);
D_Hz=diag(D/(2*pi));

m=1;
ModeSelect=0;
for i=1:length(D_Hz)
    if imag(D_Hz(i))<0 
        %negative imaginary part - ignored
    elseif imag(D_Hz(i))<1
        % modes smaller than 1Hz are not included
    elseif abs(imag(D_Hz(i))-Fbase)<1
        % modes close to base frequency are not included
    %elseif abs(real(D_Hz(i)))>20
        % with a real part that is very large
    elseif imag(D_Hz(i))>500
        % larger than 500 Hz are out of consideration
    %elseif abs(imag(D_Hz(i))) < (abs(real(D_Hz(i))) * 7)
        % larger than 30% damping ratio
    else
        ModeSelect(m)=i;
        m=m+1;
    end
end

[MdMode,ResidueAll,ZmValAll,~,~,~, ~]=...
    SimplusGT.Modal.SSCal(GminSS, N_Apparatus, ApparatusType, ModeSelect, GmDSS_Cell, GsysDSS, ApparatusInputStr, ApparatusOutputStr);

%% config C-IMR: Floating bus with maximum IMR value, SG with no effect on heat map, and IBR has the critical effect.
j=1;
clear CIMR CIMR2;
for k=1:N_Bus
    if ApparatusType{k}<30 || ApparatusType{k}>=40
        CIMR2(j).device = k;
        CIMR2(j).value = log10(100);
        CIMR2(j).mode = 0;
        j=j+1;
    end
end

%% sweep the mode
for modei=1:length(ModeSelect)
    Residue = ResidueAll{modei};
    ZmVal = ZmValAll{modei};
    SigmaMag = abs(real(MdMode(ModeSelect(modei))))*2*pi; %MdMode is in the unite of Hz, so needs to be changed to rad.

    for j=1:length(CIMR2)
        k=CIMR2(j).device;
        if ApparatusType{k} ~= 100
            IMR = SigmaMag/(SimplusGT.Frobenius_norm_dq(Residue(k))*SimplusGT.Frobenius_norm_dq(ZmVal(k)));
            IMR_o=IMR;
        if IMR<0.01
            IMR = log10(0.01);
        else
            IMR = log10(IMR);
        end

        if IMR<CIMR2(j).value
            CIMR2(j).value=IMR;
            CIMR2(j).mode = MdMode(ModeSelect(modei));
            CIMR2(j).value_orig=IMR_o;
        end
        end
    end
end

if CaseStudy<=3 
    Main_K_Plot(CIMR2,1);
else
    Main_K_Plot_14(CIMR2,1,ListLine);
end
title('Small-Signal System Strength Heatmap');


%% Calculate contribution factor
% pick out the selected mode and residue
[index,mode_pick] = searchmode_r(Modesel_,MdMode(ModeSelect));
Residue_pick = ResidueAll{index};
Zapp_pick = ZmValAll{index};
% calculate Zapp at the the selected mode
contr=zeros(Num_IBRbus,1);
realcontr_pu=zeros(Num_IBRbus,1);
imagcontr_pu=zeros(Num_IBRbus,1);
contrRealSum=0;  contrImagSum=0;
for i=1:Num_IBRbus
    contr(i) = Residue_pick(List_IBRbus(i)).dd*Zapp_pick(List_IBRbus(i)).dd + Residue_pick(List_IBRbus(i)).dq*Zapp_pick(List_IBRbus(i)).qd +...
               Residue_pick(List_IBRbus(i)).qd*Zapp_pick(List_IBRbus(i)).dq + Residue_pick(List_IBRbus(i)).qq*Zapp_pick(List_IBRbus(i)).qq;
    contrRealSum = contrRealSum + abs(real(contr(i)));
    contrImagSum = contrImagSum + abs(imag(contr(i)));
end
%normalize
realcontr_pu=real(contr)/contrRealSum;
imagcontr_pu=imag(contr)/contrImagSum;

barwidth=0.7;
figure_n2=figure_n+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.07]);
b1=bar(realcontr_pu,barwidth);
b1.EdgeColor = 'none';
axis([0 Num_IBRbus+1 -1 1]);
set(gca,'XTick',1:length(List_IBRbus));
set(gca,'XTicklabel',{' '});
set(gca,'position',[0.15,0.25,0.8,0.6]);
set(gca,'fontsize',9,'fontname','Times New Roman');

figure_n2=figure_n2+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.07]);
b2=bar(imagcontr_pu,barwidth);
b2.EdgeColor = 'none';
axis([0 Num_IBRbus+1 -1 1]);
for i=1:length(List_IBRbus)
    List_IBRbus_plot{i}=num2str(List_IBRbus(i));
end
set(gca,'XTick',1:length(List_IBRbus));
set(gca,'XTicklabel',List_IBRbus_plot);
set(gca,'position',[0.15,0.25,0.8,0.6]);
% set(gca,'fontsize',Fsize,'fontname','Times New Roman');
set(gca,'fontsize',9,'fontname','Times New Roman');


%% sensitivity analysis
pole_sel = mode_pick;  %in Hz
AppSelect = List_IBRbus; 
MdLayer3 = SimplusGT.Modal.MdLayer3(Residue_pick,Zapp_pick,mode_pick,ApparatusType,...
                AppSelect,Para,ApparatusPowerFlow,Ts,ApparatusBus,ListBus);

figure_n2=figure_n2+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.2,0.17]);
hold on
clear Imp_Para_sens Mode_Para_sens
for parasel=1:length(MdLayer3(Appselc).Result)
    Imp_Para_sens(parasel)=MdLayer3(Appselc).Result(parasel).DeltaZ;
    Mode_Para_sens(parasel)=MdLayer3(Appselc).Result(parasel).DLambdaRho_pu_Hz/real(mode_pick);
    plot([0,real(Mode_Para_sens(parasel))],[0,imag(Mode_Para_sens(parasel))]);
end
hold off
set(gca,'position',[0.15,0.2,0.7,0.7]);
axis([-12 12 -7 5])
set(gca,'XTick',-10:2:10);
set(gca,'YTick',-10:2:10);
set(gca,'fontsize',9,'fontname','Times New Roman');
return


function [index,mode_pick]=searchmode_r(mode,modelist)
distances = abs(mode - modelist);
[~, index] = min(distances(:));
mode_pick = modelist(index);
end