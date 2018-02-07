function hd = saveCurrentSetup(hd,varargin) 

% Default Parameters
    frameToSave = hd.nFrames ;
    camID = 1 ;

if ~isempty(varargin)
    frameToSave = varargin{1} ;
end

% Save Images
    if ~isempty(hd.WorkDir) && strcmp(hd.ToolBar.MainMenu.saveImages.Checked,'on')
        wd = hd.WorkDir ;
        nameImg = [wd.Path,wd.CommonName,'_',num2str(frameToSave),wd.ImagesExtension] ;
        imwrite(uint16(65535*hd.Images{frameToSave}{camID}),nameImg) ;
    end
    
% Save Data
    nameData = [wd.Path,wd.CommonName,'_InputData.mat'] ;
    InputData = hd.InputData ;
    save(nameData,'InputData') ;