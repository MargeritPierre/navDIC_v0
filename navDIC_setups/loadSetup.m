function [setup,hd] = loadSetup(hd,path)

    % SET CAMERAS AND IMAGES
        % Get Cameras folders
            camFolders = {} ;
            answer = 1 ;
            while answer~=0
                [answer] = uigetdir(path,'SELECT A CAMERA FOLDER OR CANCEL TO CONTINUE') ;
                if answer==0 ; break ; end
                camFolders{end+1} = answer ;
            end
            nCams = length(camFolders) ;
            if nCams<1
                warning('NO CAMERA DATA LOADED')
            end
        % Get associated images
            imgExt = {'png','tif','jpg','jpeg','bmp'} ;
            Images = {} ;
            Cameras = struct([]) ;
            for cam = 1:nCams    % Possible image extensions
                % Get files
                    files = dir('*.anImpossibleExtension') ;
                    for i = 1:length(imgExt)
                        f = dir([camFolders{cam},'/*.',imgExt{i}]) ;
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
                % Load images
                    wtbr = waitbar(0,'Loading Images...') ;
                    for i = 1:nImgs
                        Images{cam,i} = {imread([camFolders{cam},'/',commonName,idSTR{i},ext])} ;
                        wtbr = waitbar(i/nImgs,wtbr) ;
                    end
                    delete(wtbr) ;
                % SET THE CAMERA
                    CamName = strsplit(camFolders{cam},{'/','\'}) ;
                    Cameras(cam).Name = CamName{end} ;
                    Cameras(cam).CurrentState = 'ghost' ;
                    Cameras(cam).Adaptator = 'folder' ;
                    refImg = Images{cam,1} ;
                    Cameras(cam).VidObj.ROIPosition = [0 0 flip(size(refImg{1}))] ;
            end
            
        
    % GET INPUT DATA
        [inputFiles,inputPath] = uigetfile('*.mat','SELECT THE INPUT DATA FILE(S) OR CANCEL',path,'MultiSelect','on') ;
        if ischar(inputFiles)
            inputFiles = {inputFiles} ;
        end
        nInputs = length(inputFiles) ;
        if inputPath==0
            warning('NO INPUT DATA LOADED')
            nInputs = 0 ;
        end
        % Load files
        InputData = [] ;
        for in = 1:nInputs
            dataInFile = load([inputPath,inputFiles{in}]) ;
            for fi = fieldnames(dataInFile)
                 data = dataInFile.(fi{1}) ;
                 InputData(1:length(data),end+1) = data ;
            end
        end
            
% Put everything in the setup
    setup.Path = path ;
    setup.CommonName = commonName ;
    setup.ImagesExtension = ext ;
    
% Change the handles
    hd.Images = Images ;
    hd.InputData = InputData ;
    hd.nFrames = max(size(InputData,1),size(Images,2)) ;
    hd.CurrentFrame = 1 ;
    hd.Cameras = Cameras ;
        
    % DEFAULT TIMELINE
        hd.TimeLine = (0:hd.nFrames-1)'*[0 0 0 0 0 1] ;
    
    
    
    
    
    
    
    
    