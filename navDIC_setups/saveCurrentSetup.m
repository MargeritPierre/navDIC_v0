function hd = saveCurrentSetup(hd,varargin) 

% Default Parameters
    frameToSave = hd.nFrames ;
    camID = 1 ;

if ~isempty(varargin)
    frameToSave = varargin{1} ;
end

if ~isempty(hd.WorkDir) && strcmp(hd.ToolBar.MainMenu.saveImages.Checked,'on')
    wd = hd.WorkDir ;
    nameImg = [wd.Path,wd.CommonName,'_',num2str(frameToSave),wd.ImagesExtension] ;
    imwrite(uint16(65535*hd.Images{frameToSave}{camID}),nameImg) ;
end