classdef navDIC_CustomMacro < navDIC_AbstractMacro
%NAVDIC_CUSTOMMACRO A template for navDIC macros

properties
end

methods
    function this = navDIC_CustomMacro()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        this.Name = inputdlg('Set the name of this CustomMacro:','Macro Name',1,{this.Name}) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        disp(['Running the CustomMacro renamed ''' this.Name '''']) ;
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

