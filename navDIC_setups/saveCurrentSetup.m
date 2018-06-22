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
                if isempty(strfind('/',wd.Path))
                    folderName = [wd.Path,camsFolderName,'/',hd.Cameras(camID).Name] ;
                else
                    folderName = [wd.Path,camsFolderName,'\',hd.Cameras(camID).Name] ;
                end
                
            % Is the folder exists ?
                if ~exist(folderName,'dir') ; mkdir(folderName) ; end
            % Save the image
                if isempty(wd.ImagesExtension)
                    wd.ImagesExtension = '.tif' ;
                end
                if isempty(strfind('/',wd.Path))
                    nameImg = [folderName,'/',wd.CommonName,'_',num2str(frameToSave),wd.ImagesExtension] ;
                else
                    nameImg = [folderName,'\',wd.CommonName,'_',num2str(frameToSave),wd.ImagesExtension] ;
                end
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
                    if isempty(strfind('/',wd.Path))
                        nameData = [folderName,'/',inName,'.mat'] ;
                    else
                        nameData = [folderName,'\',inName,'.mat'] ;
                    end
                    eval([inName,' = hd.InputData(:,inID) ;']) ;
                    save(nameData,inName) ;
            end
    end
        
    
