classdef navDIC_FreeRAM < navDIC_AbstractMacro
%NAVDIC_FreeRAM Free the RAM of images

properties
    KeepLastFrames = 0 % number of frames to keep
end

methods
    function this = navDIC_FreeRAM()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        this.Name = inputdlg('Set the name of this FreeRAM:','Macro Name',1,{this.Name}) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        framesToKeep = hd.CurrentFrame + (-this.KeepLastFrames:0) ;
        for cam = 1:numel(hd.Cameras)
            emptyFrames = setdiff(1:numel(hd.Images{cam}),framesToKeep) ;
            [hd.Images{cam}{emptyFrames}] = deal([]) ;
        end
    end

    function hd = onFrameChange(this,hd)
    % Executed when the navDIC current frame changes (slider motion)
        hd = run(this,hd) ; 
    end
end

end

