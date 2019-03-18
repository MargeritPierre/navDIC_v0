function [obj,hd] = navDIC_fftdisp(obj,hd)

%disp('fftdisp')

% Params
    CorrSize = obj.corrSize ; % 40 ;

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
                imgRef = hd.Images{obj.RefFrame}{camID} ;%obj.refImgs{1} 
                PtsMov = obj.MovingPoints(:,:,frame-obj.RefFrame) ; % round() ?
                imgMov = hd.Images{frame}{camID} ;
                
            case 'rel'
                PtsRef = obj.MovingPoints(:,:,frame-1) ;
                imgRef = hd.Images{frame-1}{camID} ;
                PtsMov = obj.MovingPoints(:,:,frame-1) ;
                imgMov = hd.Images{frame}{camID} ;
        end
        if iscell(imgMov)
            imgMov = imgMov{1} ;
        end
        if iscell(imgRef)
            imgRef = imgRef{1} ;
        end
        valid = ~any(isnan(PtsMov),2) ;
        obj.MovingPoints(:,:,frame) = obj.Points*NaN ;
        obj.MovingPoints(valid,:,frame) = fftDispMethod(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize) ;
    end
    
% Compute Displacements
    obj.Displacements(:,:,frame) = obj.MovingPoints(:,:,frame)-obj.Points ;
    
% Display Object
    %disp(obj)