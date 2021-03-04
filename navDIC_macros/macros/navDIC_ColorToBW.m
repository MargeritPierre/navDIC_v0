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
    % Function executed when a new frame is added to navDIC
        for id = this.camIDs(:)'
            hd.Images{id}{hd.CurrentFrame} = sum(hd.Images{id}{hd.CurrentFrame}*(1/size(hd.Images{id}{hd.CurrentFrame},3)),3,'native') ;
        end
    end
end

end

