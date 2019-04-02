function out = prepareAnimation(hd)

    disp '------ ANIMATION -------'
    
    % Initialize
        out.Valid = false ;
    
    % Select the preview to be animated
        [IDs,valid] = selectPreviews(hd,'single') ;
        if ~valid ; warning('NO VALID PREVIEW TO EXPORT') ; return ; end
        
    % Retrieve the preview
        prev = hd.Previews{IDs} ;
        if ~prev.isValid ; warning('THE SELECTED PREVIEW IS NOT VALID') ; return ; end
        
    % Set the preview as current figure
        figure(prev.fig)
        
    % Ask for the file to be saved
        defaultName = [hd.WorkDir.Path,'\',prev.SeedName,'_',char(regexprep(string(datetime),{' ','-',':'},'_')),'.avi'] ;
        [file,path] = uiputfile('*.avi','WHERE TO SAVE THE ANIMATION ?',defaultName) ;
        if file==0 ; warning('THE SELECTED FILE IS NOT VALID') ; return ; end
        VideoFile = [path file] ;
        
    % Other parameters (has to be putted in an interface later...)
        params.VideoQuality = 75 ;
        params.FramesRecorded = [hd.nFrames:-8:1] ; %[1:1:hd.nFrames] ;
    
    % Create the Video Writer
        writerObj = VideoWriter(VideoFile) ;
        writerObj.Quality = params.VideoQuality ;
        
    % Output structure
        out.Valid = true ;
        out.prevID = IDs ;
        out.fig = prev.fig ;
        out.writerObj = writerObj ;
        out.params = params ;
            
end