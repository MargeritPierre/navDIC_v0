%% BUILD THE REF IMAGE

%% LOAD & PREPARE THE MACRO 
global hd

macro = hd.Macros(1) ;
cam = macro.Seed.CamIDs ;

%macro.setupDIC ; % re-compute all DIC maps etc.

%% BUILD REFERENCE IMAGES

RF = macro.RefFrame ;
if RF==0 ; RF = macro.Seed.RefFrame ; end

refImg = hd.Images{macro.Seed.CamIDs}(RF) ;
refImg = macro.processImgs(refImg) ;

X = macro.Seed.Points ;

IMG = cell(hd.nFrames,1) ;

wtbr = waitbar(0,'Image Difference...') ;
for fr = 1:hd.nFrames
    x = macro.Seed.MovingPoints(:,:,fr) ;
    if all(isnan(x(:))) ; break ; end
    IMG(fr) = hd.Images{macro.Seed.CamIDs}(fr) ;
    IMG(fr) = macro.processImgs(IMG(fr)) ;
    IMG{fr} = macro.transformImage(X,x,IMG{fr}) ;
    IMG{fr} = abs(IMG{fr}-refImg{1}) ;
    IMG{fr}(~macro.ROI) = 0 ;
    wtbr = waitbar(fr/hd.nFrames,wtbr) ;
end
delete(wtbr) ;
    
%% PUSH TO NAVDIC
    newCamIdx = numel(hd.Cameras)+1 ;
    newCam = hd.Cameras(cam) ;
    newCam.Name = ['DIFF | ' hd.Cameras(cam).Name] ;
    hd.Images{newCamIdx} = IMG ;
    hd.Cameras(newCamIdx) = newCam ;
    
    
    
    
    