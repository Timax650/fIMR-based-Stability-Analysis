function V=FFT_analyze(data,frequency,ts)
%FFT_ANALYZE 对扫频仿真结果进行傅里叶变换和阻抗计算，保证选取的时长是1s
%% 读取数据
Ia=data;

%% FFT并绘图
%FFT分解
Y = fft(Ia);   
L = max(size(Ia));             % Length of signal
V=Y(ceil(frequency*ts)+1,:)*2/L;
% fprintf(['Output frequency: ' num2str(ceil(frequency*ts))])




 
% L = max(size(data));  
% 
% 
% % 参数设置
% f_dist = frequency;      
% N = ceil(L/ts / f_dist);        
% 
% 
% % --- 递归 DFT 实现 ---
% % 分子系数: (2/N) * [1, 0, 0, ..., -1] (长度为 N+1)
% b = (2/N) * [1, zeros(1, N-1), -1];
% % 分母系数: [1, -exp(j*2*pi*k/N)]
% a = [1, -exp(1i * 2*pi/N)];
% dc = [1, -exp(0)];
% 
% % 滤波计算 (得到每一时刻的复数相量)
% X = filter(b, a, data);
% X_DC = filter(b, dc, data)/2;
% V = X(end,:).';
% 
% 
% % 绘图验证
% t = linspace(0,ts,L);
% subplot(2,1,1);
% amp = abs(X(:,2));
% phase = angle(X(:,2));
% plot(t, data(:,2), 'k--', t, amp, 'r', t, data(:,2)-X_DC(end,2), 'b--', 'LineWidth', 1.5);
% title('原始信号与提取幅值'); legend('含暂态信号', '递归DFT提取幅值');
% subplot(2,1,2);
% plot(t, rad2deg(phase));
% title('提取相位 (Degree)');

end

