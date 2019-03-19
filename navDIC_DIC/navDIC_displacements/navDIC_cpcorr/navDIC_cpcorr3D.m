function [obj,hd] = navDIC_cpcorr3D(obj,hd)

disp('cpcorr')

% Params
    CorrSize = obj.corrSize ; % 40 ;

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    nbCam = length(camID) ;
    %refFrame = obj.RefFrame ; % needed to know which picture is the actual start of the test

    
% ELSE, Compute DIC
    for i = 1:nbCam
        % Blah blah
        if frame<=obj.RefFrame
            obj.MovingPoints(:,:,1:frame,i) = repmat(obj.Points(:,:,i),[1 1 frame 1]) ;
            obj.Displacements(:,:,1:frame,i) = zeros([size(obj.Points,1) size(obj.Points,2) frame 1]) ;
        end

        if 0 && hd.CurrentFrame>1
            obj.RefFrame = frame;
            obj.MovingPoints(:,:,1:frame,i) = repmat(obj.Points(:,:,1),[1 1 frame 1]) ;
            obj.Displacements(:,:,1:frame,i) = zeros([size(obj.Points,1) size(obj.Points,2) frame nbCam]) ;
        end
        if frame-obj.RefFrame >= 1
            switch obj.displMode
                case 'abs'
                    PtsRef = obj.Points(:,:,i) ;
                    if iscell(hd.Images{obj.RefFrame}{i}) 
                        imgRef = hd.Images{obj.RefFrame}{i}{1} ;%obj.refImgs{1} ;
                    else
                        imgRef = hd.Images{obj.RefFrame}{i} ;
                    end
                    PtsMov = round(obj.MovingPoints(:,:,frame-obj.RefFrame,i)) ;
                    if iscell(hd.Images{frame}{i}) 
                        imgMov = hd.Images{frame}{i}{1} ;%obj.refImgs{1} ;
                    else
                        imgMov = hd.Images{frame}{i} ;
                    end
                case 'rel'
                    PtsRef = obj.MovingPoints(:,:,frame-1,i) ;
                    if iscell(hd.Images{frame-1}{i}) 
                        imgRef = hd.Images{frame-1}{i}{1} ;%obj.refImgs{1} ;
                    else
                        imgRef = hd.Images{frame-1}{i} ;
                    end
                    PtsMov = obj.MovingPoints(:,:,frame-1,i) ;
                    if iscell(hd.Images{frame}{i}) 
                        imgMov = hd.Images{frame}{i}{1} ;%obj.refImgs{1} ;
                    else
                        imgMov = hd.Images{frame}{i} ;
                    end
            end
            valid = ~any(isnan(PtsMov),2) ;
            obj.MovingPoints(:,:,frame,i) = obj.Points(:,:,i)*NaN ;
            obj.MovingPoints(valid,:,frame,i) = my_cpcorr(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize(i)) ;
        end
        % Compute Displacements
        obj.Displacements(:,:,frame,i) = obj.MovingPoints(:,:,frame,i)-obj.Points(:,:,i) ;
    end

    
% Display Object
    disp(obj)