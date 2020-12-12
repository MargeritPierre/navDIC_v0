%% ADD A LISTENER TO THE navDIC Frame Slider
global hd

if exist('lst','var') ; delete(lst) ; end

slider = hd.ToolBar.frameSlider ;

callback = @()set(txt,'String',['t = ' num2str(t(round(slider.Value)),'%.1f') ' sec']) ;

lst = addlistener(slider ...
            ,'Value','PostSet' ...
            , @(src,evt)callback() ...
            ) ;
            