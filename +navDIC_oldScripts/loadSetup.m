function [setup,hd] = loadSetup(hd,path)

    % HARD-CODED PARAMETERS (for now...)
        convertWB = true ; % Convert color images to White & Black
        normalizeLocal = false ; % Normalize each frames
        normalizeGlobal = true ; % Normalize Frame Sets (see below)

    % SET CAMERAS AND IMAGES
        disp('------ IMAGES LOADING ------')
        % Get Cameras folders
            camFolders = {} ;
            answer = 1 ;
            disp('  SELECT CAMERA FOLDERS')
            while answer~=0
                [answer] = uigetdir(path,'SELECT A CAMERA FOLDER OR CANCEL TO CONTINUE') ;
                if answer==0 ; break ; end
                camFolders{end+1} = answer ;
                disp(['    Camera ',num2str(length(camFolders)),' : ',answer])
            end
            nCams = length(camFolders) ;
            if nCams<1
                warning('NO CAMERA DATA LOADED')
            end
        % Get associated images
            imgExt = {'png','tif','tiff','jpg','jpeg','bmp','raw','nef'} ;
            Images = {} ;
            Cameras = struct([]) ;
            disp('  LOAD IMAGES')
            for cam = 1:nCams    % Possible image extensions
                % Get files
                    disp(['    Camera ',num2str(cam)])
                    files = dir('*.anImpossibleExtension') ;
                    for i = 1:length(imgExt)
                        f = dir([camFolders{cam},'/*.',imgExt{i}]) ;
                        files(end+(1:length(f))) = f ;
                    end
                % Keep file names only
                    fileNames = {files.name} ;
                    if isempty(fileNames)
                        warning(['No Valid Image Files Found in ',camFolders{cam}])
                    end
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
                    disp(['      Type: ',ext])
                % Get image ids
                    idSTR = {} ;
                    idNUM = [] ;
                    for i = 1:length(fileNames)
                        idSTR{i} = fileNames{i}(length(commonName)+1:end-length(ext)) ;
                        if ~isempty(str2double(idSTR{i}))
                            idNUM(i) = str2double(idSTR{i}) ;
                        else
                            idNUM(i) = NaN ;
                        end
                    end
                    [idNUM,ind] = sort(idNUM(~isnan(idNUM))) ;
                    idSTR = idSTR(ind(~isnan(idNUM))) ;
                    nImgs = length(idSTR) ;
                    disp(['      Frames: ',num2str(nImgs)])
                % Initialize
                    % Load function
                        loadImage = @(id) imread([camFolders{cam},'/',commonName,idSTR{id},ext]) ;
                        imData = loadImage(1) ;
                    % Get Infos
                        [nI,nJ,nColors] = size(imData) ;
                        dataType = class(imData) ;
                        disp(['      Class: ',dataType])
                    % Process function
                        processImage = {} ;
                        if convertWB % To Black & White
                            processImage{end+1} = @(img) sum(img/nColors,3,'native') ;
                            nColors = 1 ;
                        end
                        if normalizeLocal % Normalize each frame
                            processImage{end+1} = @(img) img-min(img(:)) ;
                            processImage{end+1} = @(img) img*(double(max(getrangefromclass(img)))/double(max(img(:)))) ;
                        end
                % Load images
                    Images{cam} = zeros([nI nJ nColors nImgs],dataType) ;
                    wtbr = waitbar(0,'Loading images ...') ;
                    for i = 1:nImgs
                        % Load image
                            imData = loadImage(i) ;
                        % Apply Processes
                            for p = 1:length(processImage)
                                imData = processImage{p}(imData) ;
                            end
                        % Record
                            Images{cam}(:,:,:,i) = imData ;
                        % Waitbar
                            wtbr = waitbar(i/nImgs,wtbr,['Loading images (',num2str(i),'/',num2str(nImgs),')']) ;
                    end
                    if normalizeGlobal % Normalize the entire set
                        Images{cam} = Images{cam}-min(Images{cam}(:)) ;
                        Images{cam} = Images{cam}*(double(max(getrangefromclass(Images{cam})))/double(max(Images{cam}(:)))) ;
                    end
                    delete(wtbr) ;
                % SET THE CAMERA
                    CamName = strsplit(camFolders{cam},{'/','\'}) ;
                    Cameras(cam).Name = CamName{end} ;
                    Cameras(cam).CurrentState = 'ghost' ;
                    Cameras(cam).Adaptator = 'folder' ;
                    Cameras(cam).VidObj.ROIPosition = [0 0 nJ nI] ;
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
    hd.nFrames = max(size(InputData,1),size(Images{end},4)) ;
    hd.CurrentFrame = 1 ;
    hd.Cameras = Cameras ;
        
    % DEFAULT TIMELINE
        hd.TimeLine = (0:hd.nFrames-1)'*[0 0 0 0 0 1] ;
    
    
    
    
    
    
    
    
    