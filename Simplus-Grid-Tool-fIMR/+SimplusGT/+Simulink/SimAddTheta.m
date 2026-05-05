% This function adds public static theta into the simulink model.

% Author(s): Ruiting Xu
function [Name_Thetablc,Name_Thetalabel] = SimAddTheta(Name_Model,Name_LibFile,Size_Theta,Pos_Powergui,Wbase)
% add the theta block
Name_Thetablc = 'ideal theta';
FullName_Theta = [Name_Model '/' Name_Thetablc];
add_block([Name_LibFile '/ideal theta'],FullName_Theta);

Pos_Theta = Pos_Powergui + [200,0];
set_param(FullName_Theta,'position',[Pos_Theta,Pos_Theta+Size_Theta]);
set_param(FullName_Theta,'Wbase','Wbase_sim');

set_param(gcb,'Orientation','right');

% add the lable
Name_Thetalabel = 'theta';
FullName_Thetalabel = [Name_Model '/' Name_Thetalabel];
add_block('simulink/Signal Routing/Goto',FullName_Thetalabel);
set_param(FullName_Thetalabel,'GotoTag',Name_Thetalabel);

Pos_Thetalabel = Pos_Powergui + [260,5];
set_param(FullName_Thetalabel,'position',[Pos_Thetalabel,Pos_Thetalabel+[40,20]]);
set_param(gcb,'Orientation','right');

add_line(Name_Model,...
            {[Name_Thetablc '/1']}, {[Name_Thetalabel '/1']},'autorouting','smart');
end