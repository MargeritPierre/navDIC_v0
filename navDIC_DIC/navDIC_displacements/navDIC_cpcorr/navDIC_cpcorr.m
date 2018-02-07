function [obj,hd] = navDIC_cpcorr(obj,hd)

disp('cpcorr')

% Params
    CorrSize = 10 ;

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;

% If it is the first frame
    if frame == 1
        obj.MovingPoints = obj.Points ;
        obj.Displacements = zeros(size(obj.Points)) ;
    end
    
% ELSE, Compute DIC
    if frame > 1
        switch obj.displMode
            case 'abs'
                PtsRef = obj.Points ;
                imgRef = obj.refImgs{1} ;
                PtsMov = round(obj.MovingPoints(:,:,frame-1)) ;
                imgMov = hd.Images{frame}{camID} ;
            case 'rel'
                PtsRef = obj.MovingPoints(:,:,frame-1) ;
                imgRef = hd.Images{frame-1}{camID} ;
                PtsMov = obj.MovingPoints(:,:,frame-1) ;
                imgMov = hd.Images{frame}{camID} ;
        end
        valid = ~any(isnan(PtsMov),2) ;
        obj.MovingPoints(:,:,frame) = obj.Points*NaN ;
        obj.MovingPoints(valid,:,frame) = my_cpcorr(PtsMov(valid,:),PtsRef(valid,:),imgMov,imgRef,CorrSize) ;
    end
    
% Compute Displacements
    obj.Displacements(:,:,frame) = obj.MovingPoints(:,:,frame)-obj.Points ;
    
% Display Object
    disp(obj)