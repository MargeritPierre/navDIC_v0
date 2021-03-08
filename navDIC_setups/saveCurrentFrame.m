function hd = saveCurrentFrame(hd,varargin) 
% When the Working Directory has been set, This function is executed at
% each frame 

% Default Parameters
    frameToSave = hd.nFrames ;
    camsFolderName = 'Images' ; 
    inputsFolderName = 'Inputs' ; 

% Maybe an other frame has to be saved ?
    if ~isempty(varargin)
        frameToSave = varargin{1} ;
    end

% If no workdir has been defined, return !
    if isempty(hd.WorkDir) ; return ; end
    wd = hd.WorkDir ;
    
% Save Images
    if ~isempty(hd.Cameras) && strcmp(hd.ToolBar.MainMenu.saveImages.Checked,'on')
        for camID = 1:length(hd.Cameras)
            % A folder by camera
                FolderName = fullfile(wd.Path,camsFolderName,hd.Cameras(camID).Name);
            % Is the folder exists ?
                if ~isfolder(FolderName) ; mkdir(FolderName) ; end
            % Save the image
                nameImg = sprintf([wd.CommonName,wd.ImagesExtension],frameToSave) ;
                nameImg = fullfile(FolderName,nameImg) ;
                imwrite(hd.Images{camID}{frameToSave},nameImg) ;
        end
    end
    

    % Save Images in binary files
    if ~isempty(hd.Cameras) && strcmp(hd.ToolBar.MainMenu.saveImagesBin.Checked,'on')
        
      % Are the frame data type and image dimensions already known ?
      FolderName = fullfile(wd.Path,camsFolderName);
      FileName = fullfile(FolderName,'ImagesInfo.txt');
             
     if ~isfile(FileName)
     % Save the resolution of all cameras and the data types of the frames in a text file
             NbCams = length(hd.Cameras);
        CamId = zeros(NbCams,1, 'int8');
        CamName= cell(NbCams,1);
        ImageHeight = zeros(NbCams,1, 'uint16');
        ImageWidth = zeros(NbCams,1, 'uint16');
        FrameDataType = cell(NbCams,1);
        
        for camID = 1:length(hd.Cameras)
            CamId(camID,1) = camID;
            CamName{camID,1} = hd.Cameras(camID).Name;
            ImageSize = hd.Cameras(camID).VidObj.ROIPosition;
            ImageHeight(camID,1) = ImageSize(3);
            ImageWidth(camID,1) = ImageSize(4);
            FrameDataType{camID,1} = class(hd.Images{camID}{end});
        end
             
        T = table(CamId,CamName,ImageHeight,ImageWidth,FrameDataType);
        if ~isfolder(FolderName) ; mkdir(FolderName) ; end
        writetable(T,FileName,'Delimiter',';')
     end
                

        for camID = 1:length(hd.Cameras)
            % A folder by camera
                FolderName = fullfile(wd.Path,camsFolderName,hd.Cameras(camID).Name) ;
            % Is the folder exists ?
                if ~isfolder(FolderName) ; mkdir(FolderName) ; end
            % Save the image
                nameImg = sprintf([wd.CommonName,'.bin'],frameToSave) ;
                FileName = fullfile(FolderName,nameImg) ;
                fileID = fopen(FileName,'w');
                frameDataType = class(hd.Images{camID}{frameToSave});
                fwrite(fileID,hd.Images{camID}{frameToSave},frameDataType);
                fclose(fileID);
        end

    end

    
% Save Inputs Data
    if ~isempty(hd.DAQInputs) && strcmp(hd.ToolBar.MainMenu.saveImages.Checked,'on')
        % Folder of inputs data
            FolderName =fullfile(wd.Path,inputsFolderName);
        % Is the folder exists ?
            if ~isfolder(FolderName) ; mkdir(FolderName) ; end
        % Saving... 
            for inID = 1:length(hd.DAQInputs.Inputs)
                inName = hd.DAQInputs.Inputs(inID).DataName ;
                % Save the Data
                    nameData = fullfile(FolderName,[inName,'.mat']) ;
                    eval([inName,' = hd.InputData(:,inID) ;']) ;
                    save(nameData,inName) ;
            end
        % Saving Timeline
            nameData = fullfile(FolderName,['Time','.mat']) ;
            time = sum(bsxfun(@times,bsxfun(@minus,hd.TimeLine,hd.TimeLine(1,:)),[0 0 3600*24 3600 60 1]),2) ;
            % TODO : garder l'info du temps absolu de l'acquisition !
            save(nameData, 'time') ;
            
    end
   
    
