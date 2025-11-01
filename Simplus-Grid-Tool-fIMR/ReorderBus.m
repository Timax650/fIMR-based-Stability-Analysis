% Find pure empty bus, i.e., no apparatus and no passive load
function [OrderOld2New, ApparatusSourceType]=ReorderBus(ListBus,ApparatusType,NumBus)
IndexEbus = [];
for k = 1:NumBus
    if (ListBus(k,5)==0) && (ListBus(k,6)==0) && (ListBus(k,7)==0) && (ListBus(k,8)==0)
        IndexEbus = [IndexEbus,ListBus(k,1)];
    end
end

% Get the apparatus source type:
% Notes: Remember to include other kinds of apparatuses later
for i = 1:NumBus
    if ApparatusType{i}==0 || ApparatusType{i}==1 || ApparatusType{i}==20 || ApparatusType{i}==41 || ApparatusType{i}==42 || ApparatusType{i}==28 || ApparatusType{i}==30
      	ApparatusSourceType(i) = 1;    % Voltage node
    elseif ApparatusType{i}==10 || ApparatusType{i}==11 || ApparatusType{i}==18
    	ApparatusSourceType(i) = 2;    % Current node
    elseif ApparatusType{i}==100
        ApparatusSourceType(i) = 3;    % Floating node     
    else
     	error('Error: The apparatus type is not included');
    end
end

% Based on the device source type, we get the new order
IndexVbus = find(ApparatusSourceType == 1);
IndexIbus = find(ApparatusSourceType == 2);
IndexFbus = find(ApparatusSourceType == 3);
OrderOld2New = [IndexVbus,IndexIbus,IndexFbus]; % Convert old order to new
OrderOld2NewNoFbus = [IndexVbus,IndexIbus];
for i = 1:NumBus
    OrderNew2Old(OrderOld2New(i)) = i;            % Convert new order back to old
end
for i = 1:(NumBus-length(IndexFbus))
    OrderNew2OldNoFbus_(OrderOld2New(i)) = i;
end
CounterFbus = 0;
for i = 1:length(OrderNew2OldNoFbus_)
    if OrderNew2OldNoFbus_(i) ~= 0
        OrderNew2OldNoFbus(i-CounterFbus) = OrderNew2OldNoFbus_(i);
    else
        CounterFbus = CounterFbus+1;
    end
end

% Existance of Node
if ~isempty(IndexVbus); ExistVbus = 1; else; ExistVbus = 0; end
if ~isempty(IndexIbus); ExistIbus = 1; else; ExistIbus = 0; end
if ~isempty(IndexFbus); ExistFbus = 1; else; ExistFbus = 0; end

% Re-order device source tyoe
ApparatusSourceType = ApparatusSourceType(:,OrderOld2New);

% Re-order power flow
%V = V(OrderOld2New,:);
%I = I(OrderOld2New,:);

% Re-order nodal admittance matrix
%Ybus = Ybus(OrderOld2New,OrderOld2New);
%Ybus_ = Ybus_(OrderOld2New,OrderOld2New);

% Re-order device para
% for i = 1:NumBus
%     ApparatusTypeNew{i} = ApparatusType{OrderOld2New(i)};
%     ApparatusParaNew{i} = Para{OrderOld2New(i)};
% end

end