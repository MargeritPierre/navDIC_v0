function [obj,hd] = navDIC_ICGN(obj,hd)

%disp('icgn')

% Params
    CorrSize = obj.corrSize ; % 40 ;

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    nbCam = length(hd.Cameras) ;
    %refFrame = obj.RefFrame ; % needed to know which picture is the actual start of the test
    
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
        
        % ELSE, Compute DIC
        if frame-obj.RefFrame >= 1
            switch obj.displMode
                case 'abs'
                    PtsRef = obj.Points(:,:,i) ;
                    imgRef = hd.Images{obj.RefFrame}{i} ;
                    PtsMov = obj.MovingPoints(:,:,frame-obj.RefFrame,i) ;
                    imgMov = hd.Images{frame}{i} ;
                case 'rel'
                    PtsRef = obj.MovingPoints(:,:,frame-1,i) ;
                    imgRef = hd.Images{frame-1}{i} ;
                    PtsMov = obj.MovingPoints(:,:,frame-1,i) ;
                    imgMov = hd.Images{frame}{i} ;
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
        while iscell(imgMov)
            imgMov = imgMov{1} ;
        end
        while iscell(imgRef)
            imgRef = imgRef{1} ;
        end
        valid = ~any(isnan(PtsMov),2) ;
        obj.MovingPoints(:,:,frame) = obj.Points*NaN ;
        obj.MovingPoints(valid,:,frame) = icgnCorrMethod(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize) ;
    end
    
% Compute Displacements
    obj.Displacements(:,:,frame) = obj.MovingPoints(:,:,frame)-obj.Points ;
    
% Display Object
    %disp(obj)