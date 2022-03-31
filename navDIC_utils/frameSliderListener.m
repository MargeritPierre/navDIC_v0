%% ADD A LISTENER TO THE navDIC Frame Slider

%% Declare the listener
global hd

if exist('frlst','var') ; delete(frlst) ; end

slider = hd.ToolBar.frameSlider ;

frlst = addlistener(slider ...
            ,'Value','PostSet' ...
            , @(src,evt)disp(slider.Value) ...
            ) ;
        
%% Text string update
% first create a text 'txt'
t = sum(hd.TimeLine(:,4:6).*[3600 60 1],2) ; 
t = t-t(2) ;
frlst.Callback = @(src,evt)set(txt,'String',['t = ' num2str(t(round(slider.Value)),'%.1f') ' sec']) ;

        
%% Text string update
% first create a text 'txt'
d = abs(hd.InputData(:,1)) ;
frlst.Callback = @(src,evt)set(txt,'String',[num2str(d(round(slider.Value)),'%.0f') ' Volts']) ;

        
%% Duration string update
% first create a text 'txt'
t = datetime(hd.TimeLine) ;
dt = t-t(1) ;
frlst.Callback = @(src,evt)set(txt,'String',string(dt(round(slider.Value)))) ;

%% Frame update
ax = hd.Previews(end).AxesImg ;
ctr = findobj(ax,'type','contour') ;
if isempty(ctr) ; [~,ctr] = contour(ax,hd.Images{4}{hd.CurrentFrame}) ; end
ctr.XData = flip(1:size(hd.Images{4}{end},2))+22 ;
frlst.Callback = @(src,evt)set(ctr,'ZData',hd.Images{4}{slider.Value}) ;
ctr.LevelList = linspace(ax.CLim(1),ax.CLim(2),size(unique(ax.Colormap,'rows'),1)+1) ;
ctr.LineColor = 'k' ;
ctr.ShowText = 'on' ;

%% BAR PLOT VIEW
DataToPlot = [...
                hd.InputData(:,[6 14 15]) ...
%                 hd.InputData.('"Flow rate CD"') ...
%                 hd.InputData.('"Flow rate AD2"') ...
%                 hd.InputData.('"TCP Speed"') ...
                ] ;
ratios = [.8 7.6 95] ;
clf ;
data = table2array(DataToPlot(hd.CurrentFrame,:)) ;
h = barh(data./ratios) ;
lbl = strcat(arrayfun(@(d)num2str(d,3),data,'UniformOutput',false),{'kg/min','mL/kg','mm/s'})
txt = text(h.YData+.01,h.XData,lbl,'horizontalalignment','left') ;
set(gca,'ytick',h.XData,'yticklabel',strcat(regexprep(DataToPlot.Properties.VariableNames,'"',''),'  '))
set(gca,'xlim',[0 1.25*max(table2array(DataToPlot)./ratios,[],'all')]) ;

frlst.Callback = @(src,evt)updateBarPlot(DataToPlot(round(slider.Value),:),ratios,h,txt) ;
function updateBarPlot(data,ratios,h,txt)
    data = table2array(data) ;
    lbl = strcat(arrayfun(@(d)num2str(d,3),data,'UniformOutput',false),{'kg/min','mL/kg','mm/s'}) ;
    h.YData = data./ratios ;
    set(txt(:),{'String'},lbl') ;
    pos = [h.YData(:)+.01 h.XData(:) h.XData(:)*0] ;
    set(txt(:),{'Position'},num2cell(pos,2))
end

