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
                camName = hd.Cameras(camID).Name ;
                camName = regexprep(camName,{'\(.+\)'},{''}) ;
                folderName = [wd.Path,filesep,camsFolderName,filesep,camName] ;
            % Is the folder exists ?
                if ~exist(folderName,'dir') ; mkdir(folderName) ; end
            % Save the image
                nameImg = sprintf([wd.CommonName,wd.ImagesExtension],frameToSave) ;
                nameImg = [folderName,filesep,nameImg] ;
                imwrite(hd.Images{camID}{frameToSave},nameImg) ;
        end
    end

% Save Inputs Data
    if ~isempty(hd.DAQInputs) && strcmp(hd.ToolBar.MainMenu.saveImages.Checked,'on')
        % Folder of inputs data
            folderName = [wd.Path,inputsFolderName] ;
        % Is the folder exists ?
            if ~exist(folderName,'dir') ; mkdir(folderName) ; end
        % Saving... 
            for inID = 1:length(hd.DAQInputs.Inputs)
                inName = hd.DAQInputs.Inputs(inID).DataName ;
                % Save the Data
                    nameData = [folderName,'/',inName,'.mat'] ;
                    eval([inName,' = hd.InputData(:,inID) ;']) ;
                    save(nameData,inName) ;
            end
        % Saving Timeline
            nameData = [folderName,'/','Time','.mat'] ;
            time = sum(bsxfun(@times,bsxfun(@minus,hd.TimeLine,hd.TimeLine(1,:)),[0 0 3600*24 3600 60 1]),2) ;
            save(nameData,time(:)) ;
            
    end
        
    
