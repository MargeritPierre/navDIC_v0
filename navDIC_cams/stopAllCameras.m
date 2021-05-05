function hd = stopAllCameras(hd)

    % Is there cameras ?
        if isempty(hd.Cameras) ; return ; end
        
    % Stop all cams
        for c = 1:length(hd.Cameras)
            if ~ismember(hd.Cameras(c).CurrentState,{'connected'}) ; continue ; end
            cam = hd.Cameras(c).VidObj ;
            if isstruct(cam) ; continue ; end % when the videoinput is unavailable
            if strcmp(cam.Running,'off') ; continue ; end
            stop(cam)
            while strcmp(cam.Running,'on') ; end
        end

end