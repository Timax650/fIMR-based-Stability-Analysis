
function Main_K_Plot(CIMR2,imr_case)
%% Load map data
% DataName = 'K_68Bus_IBR_Load_Data';
% DataName = 'K_68Bus_IBR_Data';
% DataName = 'K_68Bus_IBR_17_Data';
% DataName = 'K_68Bus_IBR_17_14_Data';
% DataName = 'K_68Bus_IBR_17_14_7_Data';
DataName = 'K_68Bus_IBR_18_Data';

Data = load(DataName).SaveData;

%%
KH              = Data.KH;
YbusVI          = Data.YbusVI;
YbusVIF         = Data.YbusVIF;
GbusVI          = Data.GbusVI;
GbusVIF         = Data.GbusVIF;
YbusOrigin      = Data.YbusOrigin;
%Index_Vbus      = Data.Index_Vbus;
%Index_Ibus      = Data.Index_Ibus;
%Index_Fbus      = Data.Index_Fbus;
%Index_Ebus      = Data.Index_Ebus;
Order_Old2New   = Data.Order_Old2New;
Order_New2Old   = Data.Order_New2Old;

%%
FigNum = 0;
ColorRGB();
FigSize = [0.1 0.1 0.5 0.75]*0.5;

FigNum = FigNum + 1;
figure(FigNum)
clf;
set(gcf,'units','normalized','outerposition',FigSize);

%% Plot graph
% GraphMatrix = NormMatrixElement(YbusOrigin,'DiagFlag',0);
GraphMatrix = UnitMatrixElement(YbusOrigin,'DiagFlag',0);
  %enlarge area of bus9, bus29, bus28, bus26 for better display
% K_e1=1.3; GraphMatrix(9,29)=K_e1*GraphMatrix(9,29); GraphMatrix(29,9)=K_e1*GraphMatrix(29,9);
% K_e2=10; GraphMatrix(29,28)=K_e2*GraphMatrix(29,28); GraphMatrix(28,29)=K_e2*GraphMatrix(28,29);
% K_e3=1.3; GraphMatrix(29,26)=K_e3*GraphMatrix(29,26); GraphMatrix(26,29)=K_e3*GraphMatrix(26,29);
% K_e4=1.5; GraphMatrix(28,26)=K_e4*GraphMatrix(28,26); GraphMatrix(26,28)=K_e4*GraphMatrix(26,28);
  %
GraphData = graph(GraphMatrix,'upper');
GraphFigure = plot(GraphData); grid on; hold on;
layout(GraphFigure, 'force', 'WeightEffect', 'direct'); 
highlight(GraphFigure,GraphData,'EdgeColor',[0,0,0],'LineWidth',1.1);       % Change all edges to black by default
highlight(GraphFigure,GraphData,'NodeColor',[0,0,0]);                    	% Change all nodes to black by default
highlight(GraphFigure,GraphData,'MarkerSize',4.5);
highlight(GraphFigure,GraphData,'NodeFontSize',9);
% highlight(GraphFigure,GraphData,'NodeFontWeight','bold');

%% sort out SG-bus, IBR-bus and floating bus
N_Bus = evalin('base', 'N_Bus');
ApparatusType = evalin('base', 'ApparatusType');
k1=1;
k2=1;
k3=1;
k4=1;
for k=1:N_Bus
    if ApparatusType{k} >= 10 && ApparatusType{k} <20 % GFL
        Index_GFLbus(k1)=k;
        k1=k1+1;
    elseif ApparatusType{k} >= 20 && ApparatusType{k} <30 % GFM
        Index_GFMbus(k2)=k;
        k2=k2+1;
    elseif ApparatusType{k} >= 30 && ApparatusType{k} <40 % SG
        Index_Vbus(k3)=k;
        k3=k3+1;
    elseif ApparatusType{k}==100
        Index_Fbus(k4)=k;
        k4=k4+1;
    end
end

%% Set grid-following node
highlight(GraphFigure,Index_GFLbus,'NodeColor',[0,128,0]/255);          % Change all current node to green by default

%% Set grid-forming node
highlight(GraphFigure,Index_GFMbus,'NodeColor',[0,0,255]/255);          % Change all current node to green by default

%% Set voltage node
highlight(GraphFigure,Index_Vbus,'NodeColor',[0,0,0]);      	% Change all voltage node to black by default

%% Set floating node
highlight(GraphFigure,Index_Fbus,'NodeColor',[0.7,0.7,0.7]);   	% Change all floating node to gray by default


XData_ = GraphFigure.XData';
YData_ = GraphFigure.YData';
% XData_(2) = XData_(2)+0.02;
% YData_(2) = YData_(2)+0.02;
XData_(9) = XData_(9)-0.2;
YData_(9) = YData_(9)+0.1;
XData_(26) = XData_(26)+0.2;
YData_(26) = YData_(26)+0.02;
XData_(29) = XData_(29)-0.2;
YData_(29) = YData_(29)+0.15;
XData_(28) = XData_(28)-0.1; 
YData_(28) = YData_(28)-0.05; 

XData = XData_([CIMR2(:).device]);
YData = YData_([CIMR2(:).device]);
ZData = [CIMR2(:).value].';


%% further refine the edge color
XData = [XData.', 1.3, 1.7, -3.1, 2.5, 3.6, 3.7, 1.8, 1.1, 3, 3, 3.2, 3.15].';
YData = [YData.', -0.8, -0.4, 1.1, 3.7, 2.6, 0.4, -0.9, -2.3, 0.9, 0.7 1.3, 1.25].';
% YData(28) = YData(28)+3;
ZData = [ZData.', 2, 2, 2, 2, 2, 2, 2, 2, ZData(29),ZData(29),ZData(30),ZData(30)].';
% Plot heat map
if imr_case==1 % small-signal IMR
    PlotHeatMap(XData,YData,ZData,1,[-2,0]);
    % hold on
    % plot(XData,YData);
    % hold off
elseif imr_case==2 %large-signal IMR
    PlotHeatMap(XData,YData,ZData,1,[0,1]);
end
% Get the max
ZDataMax = max(ZData);

% Move graph to top
uistack(GraphFigure,'top');
%% Set Figure Lim
FigureMargin = 0.3;
% xmax = max(abs(XData));
% ymax = max(abs(YData));
% xymax = max(xmax,ymax);
% xlim([-xymax-FigureMargin,xymax+FigureMargin]);
% ylim([-xymax-FigureMargin,xymax+FigureMargin]);
xlim([min(GraphFigure.XData)-FigureMargin,max(GraphFigure.XData)+FigureMargin]);
ylim([min(GraphFigure.YData)-FigureMargin,max(GraphFigure.YData)+FigureMargin]);

%% Save
%print(figure(1),['Graph_' DataName '.png'],'-dpng','-r600');


h = colorbar;
% h.Label.String = "log(IMR)";
h.Label.Rotation = 0;
h.Label.VerticalAlignment = "cap";
h.Label.FontWeight = 'bold';
h.Label.Position = [0.5656,2.2133,0];
h.Ticks = [-2,-1.5,-1,-0.5,0];
h.TickLabels = {' '};

%colorbar('Ticks',[-2,-1,0,1,2],...
         %'TickLabels',{'-2 very weak','-1 weak','0 normal','1 strong','2 very strong'}, 'FontWeight','bold');
% title('Small-Signal System Strength Heatmap');
end
