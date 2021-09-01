classdef navDIC_RandomImage < navDIC_VirtualCamera
%NAVDIC_RandomImage A virtual camera generating random images

properties
end

methods
    function this = navDIC_RandomImage()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        hd = setup@navDIC_VirtualCamera(this,hd) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        hd = run@navDIC_VirtualCamera(this,hd) ;
        if isempty(this.CameraID) ; return ; end
        hd.Images{this.CameraID}{hd.CurrentFrame} = permute(rand(this.ImageSize),[2 1 3]) ;
    end

    function hd = onNewFrame(this,hd)
    % Executed when a new frame is added to navDIC
        hd = run(this,hd) ; % by default, run (for backward compatibility)
    end

    function hd = onFrameChange(this,hd)
    % Executed when the navDIC current frame changes (slider motion)
    end
end

end

