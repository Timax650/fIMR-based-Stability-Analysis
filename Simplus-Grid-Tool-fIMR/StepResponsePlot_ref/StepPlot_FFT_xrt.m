%% Detuned system
clear;
clc;
close all;

sel = 1; % 1 for detuned, 2 for tuned, 3 for tuned against.
% if sel == 1
%     load('StepR1.mat','out');
% elseif sel == 2
%     load('StepR2.mat','out');
% elseif sel == 3
%     load('StepR3.mat','out');
% end
load('Step28_2.mat','out');
Ts=1/3e4*5;
t_start = 15.8;
t_end = 16.5;
d_start = ceil(t_start/Ts)+1;
d_end = ceil(t_end/Ts)+1;
tout = out.tout;
A11_S = out.ScopeData28.signals(1).values;
color1 = [1,0,0];
color2 = [0,0.45,0.74];
color3 = [0,0.5,0];
L_width=1;
Pwidth = 0.25;
Pheight = 0.8;
PInterval = 0.033;
Fsize = 15;
P1 = [0.1, 0.1,    Pwidth,   Pheight];
P2 = [0.1+Pwidth+PInterval, 0.1,   Pwidth,   Pheight];
P3 = [0.1+(Pwidth+PInterval)*2, 0.1,   Pwidth,   Pheight];
width_L=1;%2; %change wdith of the picture.

figure(1);
% set(gcf,'position',[500,500,900*1.4*width_L,200*1.4]);
% subplot('Position',P1);
set(gcf,'unit','normalized','position',[0.2,0.2,0.2,0.2]);
set(gca,'position',[0.17,0.2,0.7,0.7]);
plot(tout(d_start:d_end), A11_S(d_start:d_end),'LineWidth',L_width, 'Color',color1);
set(gca,'XTick',t_start:0.2:20);
set(gca,'XTicklabel',{'0','0.2','0.4','0.6','0.8','1.0'});
% set(gca,'XTick',0:0.2:10);
set(gca,'YTick',0:0.02:2);
axis([t_start,t_end,0.98,1.04]) ;
% axis([t_start,t_end,0.99,1.02]) ;
%title('A11 Apparent Power (pu)')
set(gca,'fontsize',Fsize,'fontname','Times New Roman');
grid on;

%FFT analysis
%% FFT并绘图
%FFT分解
t_start_FFT = 16.1;
t_end_FFT = 17.1;
d_start_FFT = ceil(t_start_FFT/Ts)+1;
d_end_FFT = ceil(t_end_FFT/Ts);
Ia = A11_S(d_start_FFT:d_end_FFT);
Y = fft(Ia);   
L = length(Ia);             % Length of signal
P = abs(Y(1:L/2+1)*2/L);      %取幅值
P(1) = abs(Y(1)/L);
Ph = angle(Y(1:L/2+1));
% Ph = angle(Y(1:L/2+1)*exp(1i*pi/2))*180/pi;    %fft函数输出的相位是基于余弦的，相位滞后正弦90度
f = 1/Ts*(0:(L/2))/L;

figure(2)
set(gcf,'unit','normalized','position',[0.2,0.2,0.16,0.12]);
set(gca,'position',[0.15,0.25,0.7,0.5]);
plot(f(1:2000),P(1:2000),'LineWidth',L_width);
axis([0 150 0 0.6e-3]);
% axis([0 150 0 0.4e-3]);
% set(gca,'XTicklabel',{' '});
set(gca,'fontsize',Fsize,'fontname','Times New Roman');