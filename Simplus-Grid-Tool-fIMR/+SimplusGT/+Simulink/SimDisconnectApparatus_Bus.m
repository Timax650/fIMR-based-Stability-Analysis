% This function connects apparatus to bus

% Author(s): Yitong Li

function SimDisconnectApparatus_Bus(Name_Model,Name_Bus,Name_Apparatus,ListBus,ApparatusType)

% Organize data
N_Apparatus = length(ListBus);

% Add block
for i = 1:N_Apparatus

    % Get the bus index of the apparatus
    Bus = ListBus(i);

    % If the apparatus is NOT a "floating bus"
    if ApparatusType{Bus}<=90
        % For ac apparatuses
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn1'],[Name_Bus{Bus} '/Lconn1']);
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn2'],[Name_Bus{Bus} '/Lconn2']);
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn3'],[Name_Bus{Bus} '/Lconn3']);
    elseif 1000<=ApparatusType{Bus} && ApparatusType{Bus}<=1090
        % For dc apparatuses
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn1'],[Name_Bus{Bus} '/Lconn1']);
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn2'],[Name_Bus{Bus} '/Lconn2']);
    elseif 2000<=ApparatusType{Bus} && ApparatusType{Bus}<=2090
        % For interlink apparatuses
        % Connect to ac bus
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn1'],[Name_Bus{Bus(1)} '/Lconn1']);
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn2'],[Name_Bus{Bus(1)} '/Lconn2']);
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn3'],[Name_Bus{Bus(1)} '/Lconn3']);
        % Connect to dc bus
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn5'],[Name_Bus{Bus(2)} '/Lconn1']);
        delete_line(Name_Model,[Name_Apparatus{Bus} '/Lconn6'],[Name_Bus{Bus(2)} '/Lconn2']);
    end

end

end