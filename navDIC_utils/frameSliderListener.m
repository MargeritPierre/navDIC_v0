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
frlst.Callback = @()set(txt,'String',['t = ' num2str(t(round(slider.Value)),'%.1f') ' sec']) ;

%% Frame update
ax = hd.Previews(end).AxesImg ;
ctr = findobj(ax,'type','contour') ;
if isempty(ctr) ; [~,ctr] = contour(ax,hd.Images{4}{hd.CurrentFrame}) ; end
ctr.XData = flip(1:size(hd.Images{4}{end},2))+22 ;
frlst.Callback = @(src,evt)set(ctr,'ZData',hd.Images{4}{slider.Value}) ;
ctr.LevelList = linspace(ax.CLim(1),ax.CLim(2),size(unique(ax.Colormap,'rows'),1)+1) ;
ctr.LineColor = 'k' ;
ctr.ShowText = 'on' ;