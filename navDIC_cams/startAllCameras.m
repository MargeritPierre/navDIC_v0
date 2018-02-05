function hd = startAllCameras(hd)

    % Is there cameras ?
        if isempty(hd.Cameras) ; return ; end
        
    % Stop all cams
        for c = 1:length(hd.Cameras)
            cam = hd.Cameras(c).VidObj ;
            cam.FramesPerTrigger = 1 ;
            cam.TriggerRepeat = Inf ;
            if strcmp(cam.Running,'on') ; continue ; end
            start(cam)
            while strcmp(cam.Running,'off') ; end
        end
        
end