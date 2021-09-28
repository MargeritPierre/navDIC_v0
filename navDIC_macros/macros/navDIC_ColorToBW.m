classdef navDIC_ColorToBW < navDIC_AbstractMacro
%NAVDIC_ColorToBW Convert colorscale images to Black and White

properties
    camIDs = [] ;
end

methods
    function this = navDIC_ColorToBW()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        [IDs,valid] = selectCameras(hd,'multiple') ;
        if ~valid ; return ; end
        this.camIDs = IDs ;
    end

    function hd = run(this,hd)
    % Updating function
        for id = this.camIDs(:)'
            hd.Images{id}{hd.CurrentFrame} = sum(hd.Images{id}{hd.CurrentFrame}*(1/size(hd.Images{id}{hd.CurrentFrame},3)),3,'native') ;
        end
    end

    function hd = onNewFrame(this,hd)
    % Executed when a new frame is added to navDIC
        hd = run(this,hd) ; % by default, run (for backward compatibility)
    end

    function hd = onFrameChange(this,hd)
    % Executed when the navDIC current frame changes (slider motion)
        hd = run(this,hd) ; 
    end
end

end

