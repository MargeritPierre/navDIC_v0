function hd = captureCameras(hd)

    % Params
        timeOut = 5 ;

    % Is there cameras ?
        if isempty(hd.Cameras) ; return ; end
        nCams = length(hd.Cameras) ;
        
    % If one Camera was stopped
        for cam = 1:nCams
            if strcmp(hd.Cameras(cam).VidObj.Running,'off')
                hd = startAllCameras(hd) ;
                return ;
            end
        end
        
    % Get all cams data
        images = cell(1,nCams) ;
        for cam = 1:nCams
            cam = hd.Cameras(cam).VidObj ;
            t = tic ;
%            while cam.FramesAvailable == 0 && toc(t)<timeOut ; end
%             if toc(t)>timeOut 
%                 images{c} = zeros(cam.VideoResolution) ;
%                 continue ; 
%             end
            images{cam} = im2single(peekdata(cam,1)) ;
        end
        
    % Flush cameras buffers
        flushdata([hd.Cameras.VidObj],'all') ;
        
    % Save Images
        for cam = 1:nCams
            hd.Images{cam}(:,:,:,hd.nFrames) = images{cam} ;
        end

end