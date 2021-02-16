classdef (Abstract) navDIC_VirtualCamera < navDIC_AbstractMacro
%NAVDIC_VirtualCamera A template for navDIC macros implementing virtual camera

properties
    CameraID = [] ;
    ImageSize = [800 640 3]
end

methods
    function this = navDIC_VirtualCamera()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        if isempty(this.CameraID) ; this.CameraID = numel(hd.Cameras)+1 ; end
        defInputs = { ...
                        'Name' , this.Name ...
                        ; 'ID' , num2str(this.CameraID) ...
                        ; 'Image Size' , mat2str(this.ImageSize) ...
                    } ;
        out = inputdlg(...
                    defInputs(:,1) ...
                    ,'Virtual Camera Properties' ...
                    ,1 ...
                    ,defInputs(:,2) ...
                    ) ;
        this.Name = out{1} ;
        this.ImageSize = str2num(out{3}) ;
        hd = this.setCam(hd,str2double(out{2})) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        disp(['         Run ' this.Name ' virtual camera']) ;
    end
end

methods
    function hd = setCam(this,hd,id)
    % Change the camera ID
    % If no ID provided
        if nargin<3
            if isempty(this.CameraID)
                id = numel(hd.Cameras)+1 ; 
            else
                id = this.CameraID ;
            end
        end
    % Delete current camera & data
        if id~=this.CameraID
            hd.Cameras(this.CameraID).CurrentState = 'deleted' ;
            hd.Images(this.CameraID) = {} ;
        end
    % Create the virtual camera
        cam = struct('Infos',[],'Adaptor','none','VidObj',[],'Name',this.Name,'CurrentState','virtual') ;
        cam.Infos.DeviceName = this.Name ;
        cam.VidObj.VideoResolution = this.ImageSize(1:2) ;
        cam.VidObj.ROIPosition = [0 0 this.ImageSize(1:2)] ;
        cam.VidObj.Running = 'on' ;
    % Add the camera to navDIC
        if isempty(hd.Cameras)
            hd.Cameras = cam ;
        else
            hd.Cameras(id) = cam ;
        end
    % Set the ID
        this.CameraID = id ;
    end
end

end

