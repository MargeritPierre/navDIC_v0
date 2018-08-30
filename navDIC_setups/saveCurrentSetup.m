function hd = saveCurrentSetup(hd,varargin) 

% Default Parameters
    frameToSave = hd.nFrames ;
    camID = 1 ;

if ~isempty(varargin)
    frameToSave = varargin{1} ;
end
    wd = hd.WorkDir ;
    
% Save Images
    if ~isempty(hd.Cameras) && ~isempty(hd.WorkDir) && strcmp(hd.ToolBar.MainMenu.saveImages.Checked,'on')
        nameImg = [wd.Path,wd.CommonName,'_',num2str(frameToSave),wd.ImagesExtension] ;
        imwrite(uint8(255*hd.Images{frameToSave}{camID}),nameImg) ;
    end
    
% Save Data
    nameData = [wd.Path,wd.CommonName,'_InputData.mat'] ;
    InputData = hd.InputData ;
    save(nameData,'InputData') ;