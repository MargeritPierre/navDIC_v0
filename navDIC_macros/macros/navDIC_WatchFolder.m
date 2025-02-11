classdef navDIC_WatchFolder < navDIC_VirtualCamera
%NAVDIC_VirtualCamera A template for navDIC macros implementing virtual camera

properties
    WatchFolder = 'E:\Partage_MEB\Temp_Movefile_image'
    Filter = '\*.bmp' ;
    Delay = 0.2 ;
    TimeOut = 30 ;
    FileList = {} ;
end

methods
    function this = navDIC_WatchFolder()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setupUI(this,hd)
    % Setup the macro (change parameters, etc)
        if isempty(this.CameraID) 
            camID = numel(hd.Cameras)+1 ; 
        else
            camID = this.CameraID ;
        end
    % SETTING UI
        defInputs = { ...
                        'Name' , this.Name ...
                        ; 'Camera ID' , num2str(camID) ...
                        ; "Watch folder" , this.WatchFolder ...
                        ; "Filter" , this.Filter ...
                        ; "Delay" , num2str(this.Delay) ...
                        ; "TimeOut" , num2str(this.TimeOut) ...
                    } ;
        out = inputdlg(...
                    defInputs(:,1) ...
                    ,'Grab MEB Macro Properties' ...
                    ,1 ...
                    ,defInputs(:,2) ...
                    ) ;
        this.Name = out{1} ;
        camID = str2double(out{2}) ;
        this.WatchFolder = out{3} ;
        this.Filter = out{4} ;
        this.Delay = str2double(out{5}) ;
        this.TimeOut = str2double(out{6}) ;
    % MACRO SETUP
        hd = setup(this,hd,camID) ;
    end

    function hd = setup(this,hd,camID)
    % MACRO SETUP FUNCTION
        if nargin<3 ; camID = this.CameraID ; end
        hd = this.setCam(hd,camID) ;
        this.FileList = this.listFolder() ; % reset the file list
        if ~isempty(this.FileList) % load existing images if needed
            answer = questdlg(...
                        "Image files are already present in the folder."...
                        +newline+"Do you want to load them ?","LOAD IMAGES ?"...
                        ,"Yes","No","Yes") ;
            if strcmpi(answer,"Yes") ; hd = this.addImages(hd,this.FileList,1) ; end
        end
    end

    function hd = onNewFrame(this,hd)%run(this,hd)
    % Function executed when a new frame is added to navDIC
        starttime = tic ;
        isnewfile = false ;
        while ~isnewfile && toc(starttime)<this.TimeOut && strcmp(hd.ToolBar.stopBtn.Visible,'on')
        % Scan the watch folder
            filelist = this.listFolder() ;
        % Has the folder content changed ?
            if numel(filelist)~=numel(this.FileList)
                isnewfile = true ;
            end
        % Draw graphical objects
            drawnow;
        end
    % If no new file appeared on the folder
        if ~isnewfile
            warning("WATCH FOLDER: No new image !") ;
            % return ;
        end
    % Delay the read to let the time to the file to we written
        pause(this.Delay) ;
    % Else, get the new file identifier
        newFiles = setdiff(filelist,this.FileList) ;
        warning("New image detected ! : "+string(newFiles))
    % Load the new image(s)
        hd = this.addImages(hd,newFiles,hd.CurrentFrame) ;
    % Update the file list
        this.FileList = filelist ;
    end
end

methods
    function filelist = listFolder(this)
    % Return the list of compatible files
        files = dir([this.WatchFolder this.Filter]) ;
        filelist = {files.name} ;
    end

    function hd = addImages(this,hd,filelist,idx0)
    % Load a collection of images in the navDIC handles
        for ii = 1:numel(filelist)
            img = imread([this.WatchFolder filesep filelist{ii}]) ;
            hd.Images{this.CameraID}{idx0+ii-1} = img ;
        end
        hd.nFrames = max(hd.nFrames,numel(hd.Images{this.CameraID})) ;
    end
end

end

