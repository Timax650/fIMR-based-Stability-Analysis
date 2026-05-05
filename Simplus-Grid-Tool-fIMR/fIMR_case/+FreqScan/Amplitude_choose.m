% function amplitude = Amplitude_choose(frequency)
% %用于生成仿真中注入扰动的幅值
% 
% amplitude_base=0.02;
% if frequency<20
%     k=0.05;
% elseif frequency>0.5 && frequency<=1
%     k=0.1;
% elseif frequency>1 && frequency<=20
%     k=0.5;
% elseif frequency>20 && frequency<=75
%     k=1;
% elseif frequency>75 && frequency<150
%     k=1;
% elseif frequency>150 && frequency<190
%     k=1;
% elseif frequency>190 && frequency<250
%     k=1;    
% elseif frequency>250 && frequency<650
%     k=1;         
% elseif frequency>650 && frequency<1300
%     k=1;
% elseif frequency>1300 
%     k=1;
% end
% amplitude=amplitude_base*k;
% 
% end
% 

function amplitude = Amplitude_choose(frequency,Y_last,i_idea)
%用于生成仿真中注入扰动的幅值

if abs(frequency)<10000
    amplitude=i_idea/abs(Y_last);
else
    amplitude=0.02;
end

if amplitude>0.05
    amplitude=0.05;
elseif amplitude<1e-4
    amplitude=1e-4;
end

% amplitude=0.02;

end