function hd = stopAllCameras(hd)

    % Is there cameras ?
        if isempty(hd.Cameras) ; return ; end
        
    % Stop all cams
        for c = 1:length(hd.Cameras)
            if strcmp(hd.Cameras(c).CurrentState,'ghost') ; continue ; end
            cam = hd.Cameras(c).VidObj ;
            if strcmp(cam.Running,'off') ; continue ; end
            stop(cam)
            while strcmp(cam.Running,'on') ; end
        end

end