classdef navDIC_ImageFolder < navDIC_VirtualCamera
%NAVDIC_ImageFolder A camera providing images from an external source

properties
    Folder char
    FileList cell 
    TimeLine (:,6)
    KeepOnlyCurrentFrame(1,1) logical = true
end

methods
    function this = navDIC_ImageFolder()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        if isempty(this.CameraID) ; this.CameraID = numel(hd.Cameras)+1 ; end
        if isempty(this.Folder) ; this.Folder = hd.WorkDir.Path ; end
    % Get the image folder
        dir = uigetdir(this.Folder,'Select the Image Folder') ;
        if dir==0 ; return ; end
        this.Folder = dir ; 
    % Camera name from the last subfolder
        sub = strsplit(this.Folder,filesep) ;
        this.Name = ['From "' sub{end} '" Folder'] ;
    % Get the file list
        this.setFileList ;
    % Deduce the image size
        if isempty(this.FileList) ; return ; end
        this.ImageSize = size(imread(this.FileList{1}),[1 2 3]) ;
    % Instanciate the cam
        hd = this.setCam(hd,this.CameraID) ;
    % Modify handles
        hd.nFrames = max(hd.nFrames,numel(this.FileList)) ;
    % Run
        hd = run(this,hd) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        if isempty(this.CameraID) ; return ; end
    % Erase the previously imported frames
        if this.KeepOnlyCurrentFrame ; hd.Images{this.CameraID} = {} ; end
    % Import the current frame
        if hd.CurrentFrame==0 ; return ; end
        if hd.CurrentFrame>numel(this.FileList) ; return ; end
        hd.Images{this.CameraID}{hd.CurrentFrame} = imread(this.FileList{hd.CurrentFrame}) ;
    end

    function hd = onNewFrame(this,hd)
    % Executed when a new frame is added to navDIC
        hd = run(this,hd) ; % by default, run (for backward compatibility)
    end

    function hd = onFrameChange(this,hd)
    % Executed when the navDIC current frame changes (slider motion)
        hd = run(this,hd) ; 
    end
    
    function setFileList(this,sorted)
    % Build the file list from compatible image file formats
        supported = imformats ;
        imgExt = [supported.ext] ; % Image files currently supported by Matlab
        files = dir('*.anImpossibleExtension') ; % Initialization
        for i = 1:length(imgExt)
            f = dir([this.Folder,'/*.',imgExt{i}]) ;
            files(end+(1:length(f))) = f ;
        end
        if isempty(files)
            warning(['No Valid Image Files Found in',path]) ; 
            return ;
        end
        this.FileList = strcat([this.Folder filesep],{files.name}) ;
    % Sort the list
        if nargin<2 ; sorted = true ; end
        if sorted ; this.FileList = sort(this.FileList) ; end
    end
    
    function setTimeLine(this)
    % Set the frame time line from files modification date
        this.TimeLine = NaN(numel(this.FileList),6) ;
        for fr = 1:numel(this.FileList)
            info = imfinfo(this.FileList{fr}) ;
            if numel(info)>1  ; info = info(1) ; end
            t = [] ;
            if isempty(t) ; try ; t = datetime(info.DateTime,'InputFormat','yyyy:MM:dd HH:mm:ss') ; end ; end
            if isempty(t) ; try ; t = datetime(info.FileModDate,'InputFormat','dd-MMM-yyyy HH:mm:ss') ; end ; end
            if isempty(t) ; try ; t = datetime(info.FileModDate,'InputFormat','dd-MMMM-yyyy HH:mm:ss','Locale','fr_FR') ; end ; end
            if isempty(t) ; try ; t = datetime(regexprep(info.FileModDate,'.\.',''),'InputFormat','dd-MMM-yyyy HH:mm:ss') ; end ; end
            if isempty(t) ; continue ; end
            this.TimeLine(fr,:) = [t.Year t.Month t.Day t.Hour t.Minute t.Second] ;
        end
    end
end

end

