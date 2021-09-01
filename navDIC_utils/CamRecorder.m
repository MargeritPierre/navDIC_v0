classdef CamRecorder < handle
%CAMRECORDER Summary of this class goes here
%   Detailed explanation goes here

properties
    VideoSource
    WriterObj
    Preview
end

methods
    function this = CamRecorder()
    % Class Constructor
        global hd
    % Select the camera
        [id,valid] = selectCameras(hd,'single') ;
        if ~valid ; warning('NO VALID PREVIEW TO EXPORT') ; return ; end
        this.VideoSource = hd.Cameras(id).VidObj ;
    % Select the video file
        defaultName = [hd.WorkDir.Path,'\','Video','_',char(regexprep(string(datetime),{' ','-',':'},'_')),'.mp4'];%,'.avi'] ;
        [file,path] = uiputfile('.mp4','WHERE TO SAVE THE ANIMATION ?',defaultName) ;
        if file==0 ; warning('THE SELECTED FILE IS NOT VALID') ; return ; end
        VideoFile = [path file] ;
    % Init the Writer object
        this.WriterObj = VideoWriter(VideoFile,'MPEG-4') ;
        this.WriterObj.Quality = 100 ;
        open(this.WriterObj) ;
    % Set the preview
        this.Preview = preview(this.VideoSource) ;
        setappdata(this.Preview,'UpdatePreviewWindowFcn',@(obj,event,hImage)updatePreview(this,obj,event,hImage)) ;
    end
    
    function updatePreview(this,src,evt,hImg)
    % On preview update
        evt
        disp('updated');
        %writeVideo(this.WriterObj,...)
    end

    function delete(this)
    % Class destructor
        if ~isempty(this.WriterObj) ; close(this.WriterObj) ; end
        if ~isempty(this.Preview) ; close(this.Preview.Parent) ; end
    end
end
end

