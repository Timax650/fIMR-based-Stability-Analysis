% Plot the calculated or measured Ysys and fIMR
freq_plot=zeros(1,length(imp_meas));
for k=1:length(ListPerturb)
    ZA_plot{k}=zeros(2,2,length(imp_meas));
    Ysys_plot{k}=zeros(2,2,length(imp_meas));
    fIMRv_plot{k}=zeros(1,1,length(imp_meas));
    ZA_theo{k}=zeros(2,2,length(imp_meas));
    Ysys_theo{k}=zeros(2,2,length(imp_meas));
    fIMRv_theo{k}=zeros(1,1,length(imp_meas));
    err_fIMR{k}=zeros(1,1,length(imp_meas));
end
for n=1:length(imp_meas)
    % for n=1:127
    freq_plot(n)=imp_meas(n).frequency;
    nc = searchmode([freq_plot(n)],OmegaPN/2/pi);
    for k=1:length(ListPerturb)
        ZA_plot{k}(:,:,n)=-imp_meas(n).ZA(:,:,k);
        Ysys_plot{k}(:,:,n)=imp_meas(n).Ysys(:,:,k);
        fIMRv_plot{k}(:,:,n)=1/norm(ZA_plot{k}(:,:,n))/norm(Ysys_plot{k}(:,:,n));
        ZA_theo{k}(:,:,n)=Zapp{k}(:,:,nc);
        Ysys_theo{k}(:,:,n)=Ysys{k}(:,:,nc);
        fIMRv_theo{k}(:,:,n)=M{k}(:,:,nc);
        err_fIMR{k}(:,:,n)=abs((fIMRv_plot{k}(:,:,n)-fIMRv_theo{k}(:,:,n))/fIMRv_plot{k}(:,:,n));
    end
end

figure_n4=2001; figure(figure_n4); clf
for k=1:Num_IBRbus
    SimplusGT.plot_c(fIMRv_plot{k}(1,1,:),freq_plot,'Marker','x','MarkerSize',4,'LineStyle','none','PhaseOn',0,'PhaseShift',0,'LineWidth',0.5);
    % SimplusGT.plot_c(ZA_plot{k}(1,1,:),freq_plot,'Marker','x','MarkerSize',4,'LineStyle','none','PhaseOn',1,'PhaseShift',0,'LineWidth',0.5);
    % SimplusGT.plot_c(Ysys_plot{k}(1,1,:),freq_plot,'Marker','x','MarkerSize',4,'LineStyle','none','PhaseOn',1,'PhaseShift',0,'LineWidth',0.5);
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

function f=freq_gen(fstart,fend,fdev)
log_f=log10(fstart):fdev:log10(fend);
f=zeros(1,length(log_f));
for k=1:length(log_f)
    % if log_f(k)>1
    %     f(k)=floor(10^(log_f(k)));
    % else
    %     f(k)=floor(10*10^(log_f(k)))/10;
    % end
    f(k)=floor(10*10^(log_f(k)))/10;
end
end