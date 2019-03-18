function [obj,hd] = navDIC_cpcorr3D(obj,hd)

disp('cpcorr')

% Params
    CorrSize = obj.corrSize ; % 40 ;

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    nbCam = length(camID) ;
    %refFrame = obj.RefFrame ; % needed to know which picture is the actual start of the test

% Blah blah
    if frame<=obj.RefFrame
        obj.MovingPoints = repmat(obj.Points(:,:,1),[1 1 frame nbCam]) ;
        obj.Displacements = zeros([size(obj.Points,1) size(obj.Points,2) frame nbCam]) ;
    end
    
    if 0 && hd.CurrentFrame>1
        obj.RefFrame = frame;
        obj.MovingPoints = repmat(obj.Points(:,:,1),[1 1 frame nbCam]) ;
        obj.Displacements = zeros([size(obj.Points,1) size(obj.Points,2) frame nbCam]) ;
    end
    
% ELSE, Compute DIC
    for i = 1:nbCam
        if frame-obj.RefFrame >= 1
            switch obj.displMode
                case 'abs'
                    PtsRef = obj.Points(:,:,i) ;
                    imgRef = hd.Images{obj.RefFrame}{i} ;%obj.refImgs{1} ;
                    PtsMov = round(obj.MovingPoints(:,:,frame-obj.RefFrame,i)) ;
                    imgMov = hd.Images{frame}{i} ;
                case 'rel'
                    PtsRef = obj.MovingPoints(:,:,frame-1,i) ;
                    imgRef = hd.Images{frame-1}{i} ;
                    PtsMov = obj.MovingPoints(:,:,frame-1,i) ;
                    imgMov = hd.Images{frame}{i} ;
            end
            valid = ~any(isnan(PtsMov),2) ;
            obj.MovingPoints(:,:,frame,i) = obj.Points*NaN ;
            obj.MovingPoints(valid,:,frame,i) = my_cpcorr(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize(i)) ;
        end
        % Compute Displacements
        obj.Displacements(:,:,frame,i) = obj.MovingPoints(:,:,frame,i)-obj.Points(:,:,i) ;
    end

    
% Display Object
    disp(obj)