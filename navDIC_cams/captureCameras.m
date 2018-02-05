function hd = captureCameras(hd)

    % Is there cameras ?
        if isempty(hd.Cameras) ; return ; end
        nCams = length(hd.Cameras) ;
        
    % If one Camera was stopped
        for c = 1:nCams
            if strcmp(hd.Cameras(c).VidObj.Running,'off')
                hd = startAllCameras(hd) ;
                return ;
            end
        end
        
    % Get all cams data
        images = cell(1,nCams) ;
        for c = 1:nCams
            cam = hd.Cameras(c).VidObj ;
            images{c} = im2single(peekdata(cam,1)) ;
        end
        
    % Flush cameras buffers
        flushdata([hd.Cameras.VidObj],'all') ;
        
    % Save Images
        hd.Images{hd.nFrames} = images ;

end