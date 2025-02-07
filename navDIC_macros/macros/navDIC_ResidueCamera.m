classdef navDIC_ResidueCamera < navDIC_VirtualCamera
%NAVDIC_VirtualCamera A template for navDIC macros implementing virtual camera

properties
    DICMacro navDIC_DICMacro = navDIC_DICMacro.empty
end

methods
    function this = navDIC_ResidueCamera()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        if isempty(this.CameraID) 
            camID = numel(hd.Cameras)+1 ; 
        else
            camID = this.CameraID ;
        end
        defInputs = { ...
                        'Name' , this.Name ...
                        ; "DIC Macro:"+newline+strjoin(strcat(string(num2cell(1:numel(hd.Macros))),':',{hd.Macros.Name}),newline) ...
                        , num2str(find(ismember(hd.Macros,this.DICMacro))) ...
                        ; 'Camera ID' , num2str(camID) ...
                    } ;
        out = inputdlg(...
                    defInputs(:,1) ...
                    ,'Residue Camera Properties' ...
                    ,1 ...
                    ,defInputs(:,2) ...
                    ) ;
        this.Name = out{1} ;
        macroID = str2double(out{2}) ;
        camID = str2double(out{3}) ;
        if ~any(ismember(macroID,1:numel(hd.Macros))) ; return ; end
        if ~isa(hd.Macros(macroID),'navDIC_DICMacro') ; return ; end
        this.DICMacro = hd.Macros(macroID) ;
        this.ImageSize = size(this.DICMacro.RefImgs{end}) ;
        hd = this.setCam(hd,camID) ;
    end

    function hd = run(this,hd)
    % Run the macro
        imgs = hd.Images{this.DICMacro.Seed.CamIDs}(hd.CurrentFrame) ;
        imgs = this.DICMacro.processImgs(imgs) ;
        X0 = this.DICMacro.Seed.MovingPoints(:,:,max(hd.CurrentFrame-1,1)) ;
        X = this.DICMacro.Seed.MovingPoints(:,:,hd.CurrentFrame) ;
        imgs = this.DICMacro.transformImage(imgs,X,X0) ;
        for ii = 1:numel(this.DICMacro.RefImgs)
            diffImg = imgs{ii} - this.DICMacro.RefImgs{ii} ;
            if isprop(this.DICMacro,'ROI') ; diffImg(~this.DICMacro.ROI) = 0 ; end
            hd.Images{this.CameraID(ii)}{hd.CurrentFrame} = diffImg ;
        end
    end
end

methods
end

end

