% participation and sensitivity analysis based on the fIMR method
% run UserMain.m first

%% enter modes and apparatus of discussion
k_wc = searchmode(List_Mode_sel,freq_plot);
fc = freq_plot(k_wc);

%% config C-fIMR: Floating bus with maximum fIMR value, SG with no effect on heat map, and IBR has the critical effect.
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

clear M_sq;
for k=1:Num_IBRbus
    M_sq(k,:)=squeeze(fIMRv_plot{k});
end

%% sweep the mode
for modei=1:length(fc)
    j_a=0;
    for j=1:length(CfIMR)
        k=CfIMR(j).device;
        if any(k == ListPerturb)
            j_a=j_a+1;
            IMR = M_sq(j_a,k_wc(modei));
            IMR_o=IMR;
            if IMR<0.01
                IMR = log10(0.01);
            else
                IMR = log10(IMR);
            end

            if IMR<CfIMR(j).value
                CfIMR(j).value=IMR;   %把这个换成IMRF
                CfIMR(j).mode = fc(modei);
                CfIMR(j).value_orig=IMR_o;
            end
        end
    end
end
 %plot heatmap
if CaseStudy<=3 
    Main_K_Plot(CfIMR,1);
else
    Main_K_Plot_14(CfIMR,1,ListLine);
end
title('Small-Signal System Strength Heatmap');

%% Calculate contribution factor
Ysys_md_val=zeros(Num_IBRbus,length(fc));
Ysys_md_abs=zeros(Num_IBRbus,length(fc));
IMR_md_val=zeros(Num_IBRbus,length(fc));
contr=zeros(Num_IBRbus,length(fc));
realcontr_pu=zeros(Num_IBRbus,length(fc));
imagcontr_pu=zeros(Num_IBRbus,length(fc));
contrRealSum=zeros(1,length(fc));
contrImagSum=zeros(1,length(fc));
for i=1:Num_IBRbus
    for j=1:length(fc)
        Ysys_md_val(i,j)=Ysys_plot{i}(1,1,k_wc(j));
        Ysys_md_abs(i,j)=abs(Ysys_md_val(i,j));
        IMR_md_val(i,j)=M_sq(i,k_wc(j));
        contr(i,j) = Ysys_plot{i}(1,1,k_wc(j))*ZA_plot{i}(1,1,k_wc(j)) + Ysys_plot{i}(2,1,k_wc(j))*ZA_plot{i}(1,2,k_wc(j)) +...
                     Ysys_plot{i}(1,2,k_wc(j))*ZA_plot{i}(2,1,k_wc(j)) + Ysys_plot{i}(2,2,k_wc(j))*ZA_plot{i}(2,2,k_wc(j));  
        contrRealSum(j) = contrRealSum(j) + abs(real(contr(i,j)));
        contrImagSum(j) = contrImagSum(j) + abs(imag(contr(i,j)));
    end
end
 %normalize
for j=1:length(fc)
    realcontr_pu(:,j)=real(contr(:,j))/contrRealSum(j);
    imagcontr_pu(:,j)=imag(contr(:,j))/contrImagSum(j);
end

No_pole_plot=No_modec;
barwidth=0.7;
figure_n2=figure_n+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.07]);
b1=bar(realcontr_pu(:,No_pole_plot),barwidth);
b1.EdgeColor = 'none';
axis([0 Num_IBRbus+1 -1 1]);
set(gca,'XTick',1:length(ListPerturb));
set(gca,'XTicklabel',{' '});
set(gca,'position',[0.15,0.25,0.8,0.6]);
set(gca,'fontsize',9,'fontname','Times New Roman');

figure_n2=figure_n2+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.22,0.07]);
b2=bar(imagcontr_pu(:,No_pole_plot),barwidth);
b2.EdgeColor = 'none';
axis([0 Num_IBRbus+1 -1 1]);
for i=1:length(ListPerturb)
    List_IBRbus_plot{i}=num2str(ListPerturb(i));
end
set(gca,'XTick',1:length(List_IBRbus));
set(gca,'XTicklabel',List_IBRbus_plot);
set(gca,'position',[0.15,0.25,0.8,0.6]);
% set(gca,'fontsize',Fsize,'fontname','Times New Roman');
set(gca,'fontsize',9,'fontname','Times New Roman');


%% sensitivity analysis
pole_sel = 1j*fc(No_modec);  %in Hz
AppSelect = ListPerturb; 
for k=1:Num_IBRbus
    Ysys_sel{k} = Ysys_plot{k}(:,:,k_wc(No_modec));
end
MdLayer3 = SimplusGT.fIMR.MdLayer3f(Ysys_sel,pole_sel,ApparatusType,...
        AppSelect,Para,ApparatusPowerFlow,Ts,ApparatusBus,ListBus);  

figure_n2=figure_n2+1;
figure(figure_n2); clf;
set(gcf,'unit','normalized','position',[0.2,0.2,0.2,0.17]);
hold on
clear Imp_Para_sens Det_Para_sens
for parasel=1:length(MdLayer3(Appselc).Result)
    Imp_Para_sens(parasel)=MdLayer3(Appselc).Result(parasel).DeltaZ;
    Det_Para_sens(parasel)=MdLayer3(Appselc).Result(parasel).Ddet;
    plot([0,real(Det_Para_sens(parasel))],[0,imag(Det_Para_sens(parasel))]);
end
hold off
set(gca,'position',[0.15,0.2,0.7,0.7]);
axis([-12 12 -7 5])
% axis([-40 40 -39 1])
set(gca,'XTick',-10:2:10);
set(gca,'YTick',-10:2:10);
set(gca,'fontsize',9,'fontname','Times New Roman');


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