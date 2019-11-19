function [obj,hd] = navDIC_fftdisp(obj,hd)

%disp('fftdisp')

% Params
    CorrSize = [] ; % obj.corrSize ; % 40 ;

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    %refFrame = obj.RefFrame ; % needed to know which picture is the actual start of the test

% Blah blah
    if frame<=obj.RefFrame
%         obj.MovingPoints = repmat(obj.Points,[1 1 frame]) ;
%         obj.Displacements = zeros([size(obj.Points) frame]) ;
        obj.MovingPoints(:,:,frame) = obj.Points ;
%         obj.Displacements(:,:,frame) = 0*obj.Points ;
        obj.MovingPoints_bkp = obj.MovingPoints ;
    end
    
    if 0 && hd.CurrentFrame>1
        obj.RefFrame = frame;
        obj.MovingPoints = repmat(obj.Points,[1 1 frame]) ;
%         obj.Displacements = zeros([size(obj.Points) frame]) ;
    end
    
% ELSE, Compute DIC
    if frame-obj.RefFrame >= 1
        % Reference frame and points
            switch obj.displMode
                case 'abs'
                    PtsRef = obj.Points ;
                    imgRef = hd.Images{camID}(:,:,:,obj.RefFrame) ;
                case 'rel'
                    PtsRef = obj.MovingPoints(:,:,frame-1) ;
                    imgRef = hd.Images{camID}(:,:,:,frame-1) ;
            end
        % Current Frame and points
            PtsMov = obj.MovingPoints(:,:,frame-1) ;
            if obj.useExistingDisp && size(obj.MovingPoints,3)>=frame
                existPts = obj.MovingPoints(:,:,frame) ;
                useable = ~any(isnan(existPts),2) ;
                PtsMov(useable,:) = existPts(useable,:) ;
            end
            imgMov = hd.Images{camID}(:,:,:,frame) ;
        % Compute new positions
            valid = ~any(isnan(PtsMov),2) ;
            obj.MovingPoints(:,:,frame) = obj.Points*NaN ;
            obj.MovingPoints(valid,:,frame) = fftDispMethod(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize) ;
    end
    
% Compute Displacements
%     obj.Displacements(:,:,frame) = obj.MovingPoints(:,:,frame)-obj.Points ;
    
% Compute DataFields
    obj.computeDataFields ;
    
% Display Object
    %disp(obj)