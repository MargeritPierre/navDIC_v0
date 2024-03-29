classdef navDIC_Mono12To16 < navDIC_AbstractMacro
%NAVDIC_Mono12To16 Convert 12 bit images to 16 bit images

properties
    camID = []
end

methods
    function this = navDIC_Mono12To16()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        this.Name = inputdlg('Set the name of this 12To16:','Macro Name',1,{this.Name}) ;
        if ~isempty(hd.Cameras) ; this.camID = numel(hd.Cameras) ; end
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        if isempty(this.camID) ; return ; end
        hd.Images{this.camID}{hd.CurrentFrame} = uint16(16*double(hd.Images{this.camID}{hd.CurrentFrame})) ;
    end
end

end

