%Author: Yue Zhu, Ruiting Xu
%Introduction: in this file, parameters of an apparatus Ak is pertrubed one by
%one. The perturbation amplitude Delta_rho = 1e-5*(1+abs(rho)).
% At this toolbox's version, the apparatus's parameters will not affect the 
% results of power flow, therefore, there is no need to calculate power flow again.
% Eventhough the change of rho will affect the equibilium of Ak, it will
% not affect the equilibrium point at other points. Therefore, to calculate
% the new impedance, we only need to update the parameter, and then call the function
% SimplusGT.Toolbox.ApparatusModelCreate, and take the output GmDSS_Cell.
function Layer3Result = MdLayer3f(Ysys,Mode_Hz,ApparatusType,...
                ApparatusSelL3All,Para,ApparatusPowerFlow,Ts,ApparatusBus,ListBus)

ApparatusSelNum=length(ApparatusSelL3All);
Mode_rad = Mode_Hz*2*pi;
for ApparatusCount = 1:ApparatusSelNum
    ApparatusSelL3 = ApparatusSelL3All(ApparatusCount);
    ParamName = fieldnames(Para{ApparatusSelL3});
    ParamNum = length(ParamName);
    Ysys_ = Ysys{ApparatusCount};
    %perturb the parameters one by one.
    for k=1:ParamNum
        % original para
        [~,GmDSS_Cell_Orig,~,~,~,~,~,~,~] ...
            = SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{ApparatusSelL3},ApparatusType{ApparatusSelL3},...
            ApparatusPowerFlow{ApparatusSelL3},Para{ApparatusSelL3},Ts,ListBus);
        ZmValOrig = SimplusGT.Modal.ApparatusImpedanceCal(GmDSS_Cell_Orig, Mode_rad);

        % new para
        ParaSel = getfield(Para{ApparatusSelL3},ParamName{k}); % extract the parameter
        delta_para = 1e-5*(1e-2+abs(ParaSel));
        ParaPerturb = ParaSel + delta_para ; % add perturabation
        ParaNew = setfield(Para{ApparatusSelL3}, ParamName{k}, ParaPerturb); % update the parameter  
   
        [~,GmDSS_Cell_New,~,~,~,~,~,~,~] ...
        = SimplusGT.Toolbox.ApparatusModelCreate(ApparatusBus{ApparatusSelL3},ApparatusType{ApparatusSelL3},...
                            ApparatusPowerFlow{ApparatusSelL3},ParaNew,Ts,ListBus); 
        ZmValNew = SimplusGT.Modal.ApparatusImpedanceCal(GmDSS_Cell_New, Mode_rad);
        
        Layer3Result(ApparatusCount).Apparatus={['Apparatus',num2str(ApparatusSelL3)]};
        Layer3Result(ApparatusCount).Result(k).ParaName = ParamName(k);
        Layer3Result(ApparatusCount).Result(k).DeltaZ.dd = (ZmValNew.dd - ZmValOrig.dd)/(delta_para);
        Layer3Result(ApparatusCount).Result(k).DeltaZ.dq = (ZmValNew.dq - ZmValOrig.dq)/(delta_para);
        Layer3Result(ApparatusCount).Result(k).DeltaZ.qd = (ZmValNew.qd - ZmValOrig.qd)/(delta_para);
        Layer3Result(ApparatusCount).Result(k).DeltaZ.qq = (ZmValNew.qq - ZmValOrig.qq)/(delta_para);

        Layer3Result(ApparatusCount).Result(k).Ddet = ...
            Layer3Result(ApparatusCount).Result(k).DeltaZ.dd * Ysys_(1,1) * abs(ParaSel)...
            + Layer3Result(ApparatusCount).Result(k).DeltaZ.dq *  Ysys_(2,1) * abs(ParaSel)...
            + Layer3Result(ApparatusCount).Result(k).DeltaZ.qd *  Ysys_(1,2) * abs(ParaSel)...
            + Layer3Result(ApparatusCount).Result(k).DeltaZ.qq *  Ysys_(2,2) * abs(ParaSel);
        
    end
end
end