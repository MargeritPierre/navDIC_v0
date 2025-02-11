classdef (Abstract) navDIC_AbstractMacro < handle ...
                                            & matlab.mixin.Heterogeneous ...
                                            & matlab.mixin.Copyable
%NAVDIC_ABSTRACTMACRO Base class for navDIC macros

properties
    Name char % the macro name
    Enable logical = true
end

methods
    function this = navDIC_AbstractMacro()
    % Class constructor
        if isempty(this.Name) ; this.Name = regexprep(class(this),'navDIC_','') ; end
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setupUI(this,hd)
    % Setup the macro via its UI (change parameters, etc)
        hd = setup(this,hd) % for backward compatibility
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
    end

    function hd = run(this,hd)
    % Run the macro
    end

    function hd = onNewFrame(this,hd)
    % Executed when a new frame is added to navDIC
        hd = run(this,hd) ; % by default, run (for backward compatibility)
    end

    function hd = onFrameChange(this,hd)
    % Executed when the navDIC current frame changes (slider motion)
    end
end

methods (Sealed=true)
    function tf = eq(varargin)
        tf = eq@handle(varargin{:});
    end
end

end

