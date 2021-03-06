function hd = startAllCameras(hd)

    % Is there cameras ?
        if isempty(hd.Cameras) ; return ; end
        
    % Start all cams
        for c = 1:length(hd.Cameras)
            if ~ismember(hd.Cameras(c).CurrentState,{'connected'}) ; continue ; end
            cam = hd.Cameras(c).VidObj ;
            cam.FramesPerTrigger = 1 ;
            cam.TriggerRepeat = Inf ;
            if strcmp(cam.Running,'on') ; continue ; end
            start(cam)
            while strcmp(cam.Running,'off') ; end
        end
        
end