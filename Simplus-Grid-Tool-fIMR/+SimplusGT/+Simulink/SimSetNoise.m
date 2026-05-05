function nset = SimSetNoise(Noise_enable,Noise_raw,PowerFlow,ListPerturb,Fs)
% Initializing noise blocks in the simulink model
seed_raw = randi(10^6, 4, length(ListPerturb));

% current noise setting
for k=1:length(ListPerturb)
    nset{k}.seed = seed_raw(:,k);
    amp_cur_NI{k} = zeros(1,length(Noise_raw));
    normalized_f_NI{k} = zeros(1,length(Noise_raw));
    for n=1:length(Noise_raw)
        amp_cur_NI{k}(n) = Noise_raw(n).magnitude;
        normalized_f_NI{k}(n) = Noise_raw(n).frequency;
    end
    cur_co(k)=sqrt(PowerFlow{ListPerturb(k)}(1)^2+PowerFlow{ListPerturb(k)}(2)^2)+0.01;
    amp_cur_NI{k}=amp_cur_NI{k}*cur_co(k)/sqrt(2);
    normalized_f_NI{k} = 1 * normalized_f_NI{k} /Fs;

    % Calculate coefficients of noise blocks
    nset{k}.co_cur = fir2(30, normalized_f_NI{k}, amp_cur_NI{k});
    nset{k}.co_vol = 0.02/sqrt(2);
    nset{k}.enable = Noise_enable;
end