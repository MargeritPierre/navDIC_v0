classdef navDIC_SaveFramesToHardDrive < navDIC_AbstractMacro
%NAVDIC_CUSTOMMACRO A template for navDIC macros

properties
end

methods
    function this = navDIC_SaveFramesToHardDrive()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
        %this.Name = inputdlg('Set the name of this CustomMacro:','Macro Name',1,{this.Name}) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        disp(['Running the CustomMacro renamed ''' this.Name '''']) ;
        
    % Is there cameras ?
        if isempty(hd.Cameras) ; return ; end
        nCams = length(hd.Cameras) ;
        
    % Erase the latest image acquired and replace it by an empty cell
        for c = 1:nCams
            nFrames = length(hd.Images{c});
            if nFrames > 1
                hd.Images{c}{nFrames - 1} = {};
            end
        end
    end
end

end
