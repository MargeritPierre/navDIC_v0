function [obj,hd] = navDIC_fftdisp(obj,hd)

%disp('fftdisp')

% Params
    CorrSize = [] ; % obj.corrSize ; % 40 ;

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    %refFrame = obj.RefFrame ; % needed to know which picture is the actual start of the test

% Depending on the frame...
    if frame<=obj.RefFrame || frame < min(obj.FrameRange)
        obj.MovingPoints(:,:,frame) = obj.Points ;
    elseif frame > max(obj.FrameRange)
        obj.MovingPoints(:,:,frame) = obj.MovingPoints(:,:,max(obj.FrameRange)) ;
    elseif frame-obj.RefFrame >= 1 % Compute DIC
        % Reference frame and points
            if obj.EulerianReference
                PtsRef = obj.Points ;
                imgRef = hd.Images{camID}{frame-1} ;
            else
                switch obj.displMode
                    case 'abs'
                        PtsRef = obj.Points ;
                        imgRef = obj.refImgs{1} ; hd.Images{camID}{obj.RefFrame} ;
                    case 'rel'
                        PtsRef = obj.MovingPoints(:,:,frame-1) ;
                        imgRef = hd.Images{camID}{frame-1} ;
                end
            end
        % Current Frame and points
            PtsMov = obj.MovingPoints(:,:,frame-1) ;
            if obj.useExistingDisp && size(obj.MovingPoints,3)>=frame
                existPts = obj.MovingPoints(:,:,frame) ;
                useable = ~any(isnan(existPts),2) ;
                PtsMov(useable,:) = existPts(useable,:) ;
            end
            imgMov = hd.Images{camID}{frame} ;
        % Compute new positions
            valid = ~any(isnan(PtsMov),2) ;
            obj.MovingPoints(:,:,frame) = obj.Points*NaN ;
            obj.MovingPoints(valid,:,frame) = fftDispMethod(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize) ;
    end
    
% Compute DataFields
    obj.computeDataFields([],hd.CurrentFrame) ;
    %obj.computeDataFields ;
    
% Display Object
    %disp(obj)