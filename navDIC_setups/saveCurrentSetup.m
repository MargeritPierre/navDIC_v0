function hd = saveCurrentSetup(hd,varargin) 

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
                folderName = [wd.Path,camsFolderName,'/',hd.Cameras(camID).Name] ;
            % Is the folder exists ?
                if ~exist(folderName,'dir') ; mkdir(folderName) ; end
            % Save the image
                nameImg = [folderName,'\',wd.CommonName,'_',num2str(frameToSave),wd.ImagesExtension] ;
                imwrite(uint8(255*hd.Images{frameToSave}{camID}),nameImg) ;
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
            save(nameData,'time') ;
            
    end
        
    
