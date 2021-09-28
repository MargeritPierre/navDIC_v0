classdef navDIC_ImageProcessing < navDIC_AbstractMacro
%NAVDIC_ImageProcessing Convert colorscale images to Black and White

properties
    camIDs = [] ;
    BlackAndWhite logical = false
    Normalize logical = false
    Multiply double = 1
end

methods
    function this = navDIC_ImageProcessing()
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
    % Get options
        defInputs = { ...
                        'Black and White' , num2str(this.BlackAndWhite) ...
                        ; 'Normalize' , num2str(this.Normalize) ...
                        ; 'Multiply' , num2str(this.Multiply) ...
                    } ;
        out = inputdlg(defInputs(:,1),'Image Process Options',1,defInputs(:,2)) ;
        if isempty(out) ; return ; end
        this.BlackAndWhite = str2double(out{1}) ;
        this.Normalize = str2double(out{2}) ;
        this.Multiply = str2double(out{3}) ;
    end

    function hd = run(this,hd)
    % Updating function
        for id = this.camIDs(:)'
            img = hd.Images{id}{hd.CurrentFrame} ;
            if this.BlackAndWhite ; img = sum(img*(1/size(img,3)),3,'native') ; end
            if this.Normalize ; img = img*(max(getrangefromclass(img))/double(max(img(:)))) ; end
            if this.Multiply~=1 ; img = img*this.Multiply ; end
            hd.Images{id}{hd.CurrentFrame} = img ;
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

