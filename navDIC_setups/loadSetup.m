function out = loadSetup(path)
    
% Image files present in the directory (in the format: commonName_id.ext)
    % Possible image extensions
        imgExt = {'tif','png','bmp','jpg','jpeg'} ;
    % Get files
        files = dir('*.anImpossibleExtension') ;
        for i = 1:length(imgExt)
            f = dir([path,'/*.',imgExt{i}]) ;
            files(end+(1:length(f))) = f ;
        end
    % Keep file names only
        fileNames = {files.name} ;
    % Get the common name and extension
        str = strsplit(fileNames{1},'_') ;
        if length(str)==1
            commonName = str{1} ;
        else
            commonName = strjoin(str(1:end-1),'_') ;
        end
        [~,~,ext] = fileparts(str{end}) ;
    % Get image ids
        idSTR = {} ;
        idNUM = [] ;
        for i = 1:length(fileNames)
            idSTR{i} = fileNames{i}(length(commonName)+2:end-length(ext)) ;
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
    images = {} ;
    for i = 1:nImgs
        images{i} = imread([path,'/',commonName,'_',idSTR{i},ext]) ;
    end
    
% Load ROI
    ROI = [] ;
    try
        ROI = imread([path,'/',commonName,'_ROI',ext]) ;
        ROI = uint8(logical(ROI(:,:,1))) ;
    end
    
% Load DIC
    DIC = {} ;
    try
        DIC = load([path,'/',commonName,'_DIC.mat']) ;
    end
    
% Load Supplementary Data
    Data = [] ;
    try
        load([path,'/',commonName,'_data.mat']) ;
    end
    
% Put everything in output
    out.Path = path ;
    out.CommonName = commonName ;
    out.ImagesExtension = ext ;
    out.Images = images ;
    out.ROI = ROI ;
    out.DIC = DIC ;
    out.Data = Data ;