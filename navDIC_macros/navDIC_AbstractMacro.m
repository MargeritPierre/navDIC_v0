classdef (Abstract) navDIC_AbstractMacro < handle ...
                                            & matlab.mixin.Heterogeneous ...
                                            & matlab.mixin.Copyable
%NAVDIC_ABSTRACTMACRO Base class for navDIC macros

properties
    Name char % the macro name
end

methods
    function this = navDIC_AbstractMacro()
    % Class constructor
        if isempty(this.Name) ; this.Name = regexprep(class(this),'navDIC_','') ; end
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
    end
end

end

