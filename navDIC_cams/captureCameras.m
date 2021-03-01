function hd = captureCameras(hd)

    % Params
        timeOut = 5 ;

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
            t = tic ;
%            while cam.FramesAvailable == 0 && toc(t)<timeOut ; end
%             if toc(t)>timeOut 
%                 images{c} = zeros(cam.VideoResolution) ;
%                 continue ; 
%             end
            images{c} = peekdata(cam,1) ;
                
            if cam.VideoFormat == 'Mono12'
                images{c} = (2^16/2^12)*images{c};
            end
            %images{c} = im2single(images{c}) ;
        end
        
    % Flush cameras buffers
        flushdata([hd.Cameras.VidObj],'all') ;
        
    % Save Images
        for c = 1:nCams
            hd.Images{c}{hd.nFrames} = images{c} ;
        end

end