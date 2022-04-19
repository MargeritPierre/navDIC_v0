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
            if ~ismember(hd.Cameras(c).CurrentState,{'connected'}) ; continue ; end
            cam = hd.Cameras(c).VidObj ;
            t = tic ;
%            while cam.FramesAvailable == 0 && toc(t)<timeOut ; end
%             if toc(t)>timeOut 
%                 images{c} = zeros(cam.VideoResolution) ;
%                 continue ; 
%             end
            %images{c} = getdata(cam,1) ;
            images{c} = peekdata(cam,1) ;
            flushdata(cam) ; 
            %images{c} = im2single(images{c}) ;
        end
        
    % Flush cameras buffers
        isConnected = ismember({hd.Cameras.CurrentState},{'connected'}) ;
        if any(isConnected) ; flushdata([hd.Cameras(isConnected).VidObj],'all') ; end
        
    % Save Images
        for c = 1:nCams
            hd.Images{c}{hd.nFrames} = images{c} ;
        end

end