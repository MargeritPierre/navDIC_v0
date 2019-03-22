function [obj,hd] = navDIC_cpcorr(obj,hd)

disp('cpcorr')

% Params
try
    CorrSize = obj.corrSize ; % 40 ;
catch
    CorrSize = 6 ;
end
% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    %refFrame = obj.RefFrame ; % needed to know which picture is the actual start of the test

% Blah blah
    if frame<=obj.RefFrame
        obj.MovingPoints = repmat(obj.Points,[1 1 frame]) ;
        obj.Displacements = zeros([size(obj.Points) frame]) ;
    end
    
    if 0 && hd.CurrentFrame>1
        obj.RefFrame = frame;
        obj.MovingPoints = repmat(obj.Points,[1 1 frame]) ;
        obj.Displacements = zeros([size(obj.Points) frame]) ;
    end
    
% ELSE, Compute DIC
    if frame-obj.RefFrame >= 1
        switch obj.displMode
            case 'abs'
                PtsRef = obj.Points ;
                if iscell(hd.Images{obj.RefFrame}{camID})
                    imgRef = hd.Images{obj.RefFrame}{camID}{1} ;%obj.refImgs{1} ;
                else
                    imgRef = hd.Images{obj.RefFrame}{camID} ;
                end
                PtsMov = round(obj.MovingPoints(:,:,frame-obj.RefFrame)) ;
                if iscell(hd.Images{frame}{camID})
                    imgMov = hd.Images{frame}{camID}{1} ;
                else
                    imgMov = hd.Images{frame}{camID} ;
                end
            case 'rel'
                PtsRef = obj.MovingPoints(:,:,frame-1) ;
                if iscell(hd.Images{frame-1}{camID})
                    imgRef = hd.Images{frame-1}{camID}{1} ;%obj.refImgs{1} ;
                else
                    imgRef = hd.Images{frame-1}{camID} ;
                end
                PtsMov = obj.MovingPoints(:,:,frame-1) ;
                if iscell(hd.Images{frame}{camID})
                    imgMov = hd.Images{frame}{camID}{1} ;
                else
                    imgMov = hd.Images{frame}{camID} ;
                end
        end
        valid = ~any(isnan(PtsMov),2) ;
        obj.MovingPoints(:,:,frame) = obj.Points*NaN ;
        obj.MovingPoints(valid,:,frame) = my_cpcorr(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize) ;
    end
    
% Compute Displacements
    obj.Displacements(:,:,frame) = obj.MovingPoints(:,:,frame)-obj.Points ;
    
% Display Object
    disp(obj)