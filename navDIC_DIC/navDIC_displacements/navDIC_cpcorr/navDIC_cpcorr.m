function [obj,hd] = navDIC_cpcorr(obj,hd)

disp('cpcorr')

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
                imgRef = hd.Images{camID}{obj.RefFrame} ;%obj.refImgs{1} ;
                PtsMov = round(obj.MovingPoints(:,:,frame-obj.RefFrame)) ;
                imgMov = hd.Images{camID}{frame} ;
            case 'rel'
                PtsRef = obj.MovingPoints(:,:,frame-1) ;
                imgRef = hd.Images{camID}{frame-1} ;
                PtsMov = obj.MovingPoints(:,:,frame-1) ;
                imgMov = hd.Images{camID}{frame} ;
        end
        valid = ~any(isnan(PtsMov),2) ;
        obj.MovingPoints(:,:,frame) = obj.Points*NaN ;
        obj.MovingPoints(valid,:,frame) = my_cpcorr(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize) ;
    end
    
% Compute Displacements
    obj.Displacements(:,:,frame) = obj.MovingPoints(:,:,frame)-obj.Points ;
    
% Display Object
    disp(obj)