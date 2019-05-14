function hd = captureCameras(hd)

    % Params
        timeOut = 5 ;

    % Is there cameras ?
        if isempty(hd.Cameras); return ; end
        nCams = length(hd.Cameras) ;
        
    % If one Camera was stopped
        for c = 1:nCams
            if ~strcmp(hd.Cameras.CurrentState,'ghost') && strcmp(hd.Cameras(c).VidObj.Running,'off') 
                hd = startAllCameras(hd) ;
                return ;
            end
        end
        
        
    % Get all cams data
        images = cell(1,nCams) ;
        imgExt = {'png','tif','jpg','jpeg','bmp'} ;
        for c = 1:nCams
            if hd.Cameras(c).CurrentState == 'ghost'
                t = tic ;
                files = dir('*.anImpossibleExtension') ;
                for i = 1:length(imgExt)
                    f = dir([hd.Cameras(c).Path,'/*.',imgExt{i}]) ;
                    files(end+(1:length(f))) = f ;
                end
                % Keep file names only
                fileNames = {files.name} ;
                % Get the common name and extension
                str = strsplit(fileNames{1},'_') ;
                if length(str)==1
                    % Empty common name ?
                        if strcmp(str,fileNames{1})
                            commonName = '' ;
                        else
                            commonName = str{1} ;
                        end
                else
                    commonName = [strjoin(str(1:end-1),'_'),'_'] ;
                end
                [~,~,ext] = fileparts(str{end}) ;
                % Get image ids
                idSTR = {} ;
                idNUM = [] ;
                for i = 1:length(fileNames)
                    idSTR{i} = fileNames{i}(length(commonName)+1:end-length(ext)) ;
                    if ~isempty(str2num(idSTR{i}))
                        idNUM(i) = str2num(idSTR{i}) ;
                    else
                        idNUM(i) = NaN ;
                    end
                end
                [idNUM,ind] = sort(idNUM(~isnan(idNUM))) ;
                idSTR = idSTR(ind(~isnan(idNUM))) ;
                nImgs = length(idSTR) ;
                images{c} = {imread([hd.Cameras(c).Path,'/',commonName,idSTR{nImgs},ext])} ;
            else
                cam = hd.Cameras(c).VidObj ;
                t = tic ;
%            while cam.FramesAvailable == 0 && toc(t)<timeOut ; end
%             if toc(t)>timeOut 
%                 images{c} = zeros(cam.VideoResolution) ;
%                 continue ; 
%             end
                images{c} = im2single(peekdata(cam,1)) ;
        end
        
    % Flush cameras buffers
        if ~strcmp(hd.Cameras.CurrentState,'ghost')
            flushdata([hd.Cameras.VidObj],'all') ;
        end
    % Save Images
        hd.Images{hd.nFrames} = images ;

end