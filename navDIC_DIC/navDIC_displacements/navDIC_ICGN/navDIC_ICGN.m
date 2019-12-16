function [obj,hd] = navDIC_ICGN(obj,hd)

%disp('icgn')

% Params
    CorrSize = obj.corrSize ; % 40 ;

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camIDs = obj.CamIDs ;
    nbCam = length(camIDs) ;
    %refFrame = obj.RefFrame ; % needed to know which picture is the actual start of the test
    
    for i = 1:nbCam
        % Blah blah
        if frame<=obj.RefFrame 
            obj.MovingPoints(:,:,1:frame,i) = repmat(obj.Points(:,:,i),[1 1 frame 1]) ;
            obj.Displacements(:,:,1:frame,i) = zeros([size(obj.Points,1) size(obj.Points,2) frame 1]) ;
        end

        if 0 && hd.CurrentFrame>1
            obj.RefFrame = frame;
            obj.MovingPoints(:,:,1:frame,i) = repmat(obj.Points(:,:,i),[1 1 frame 1]) ;
            obj.Displacements(:,:,1:frame,i) = zeros([size(obj.Points,1) size(obj.Points,2) frame nbCam]) ;
        end
        
        % ELSE, Compute DIC
        if frame-obj.RefFrame >= 1
            switch obj.displMode
                case 'abs'
                    PtsRef = obj.Points(:,:,i) ;
                    imgRef = hd.Images{obj.RefFrame}{camIDs(i)} ;
                    imgMov = hd.Images{frame}{camIDs(i)} ;
                    disp(num2str(frame));
                    if isnan(sum(sum(obj.MovingPoints(:,:,frame,i))))
                        PtsMov = obj.MovingPoints(:,:,frame-1,i) ;
                    else
                        PtsMov = obj.MovingPoints(:,:,frame,i) ;
                    end
                case 'rel'
                    PtsRef = obj.MovingPoints(:,:,frame-1,i) ;
                    imgRef = hd.Images{frame-1}{camIDs(i)} ;
                    PtsMov = obj.MovingPoints(:,:,frame-1,i) ;
                    imgMov = hd.Images{frame}{camIDs(i)} ;
            end
            while iscell(imgMov)
                imgMov = imgMov{1} ;
            end
            while iscell(imgRef)
                imgRef = imgRef{1} ;
            end
            valid = ~any(isnan(PtsMov),2) ;
            obj.MovingPoints(:,:,frame,i) = obj.Points(:,:,i)*NaN ;
            if isa(imgMov,'uint8') 
                imgMov = double(imgMov)/255 ;
            end
            if isa(imgRef,'uint8') 
                imgRef = double(imgRef)/255 ;
            end
            obj.MovingPoints(valid,:,frame,i) = icgnCorrMethod(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize(i,:)) ;
        end
        % Compute Displacements
        obj.Displacements(:,:,frame,i) = obj.MovingPoints(:,:,frame,i)-obj.Points(:,:,i) ;
    end
    
% Display Object
    %disp(obj)