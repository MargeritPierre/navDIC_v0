function [obj,hd] = navDIC_deltaL(obj,hd)

disp('deltaL_L')
    
% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    PtsMov = round(obj.MovingPoints(:,:,frame)) ; % 
    Pts = obj.Points; % ref absolu dl/l0
    PtsMovRef = round(obj.MovingPoints(:,:,frame-1)); % ref relatif dl/l  
    disp_Mode = obj.displMode ;
    nPts = length(Pts(:,:,frame)) ;
    L0 = ones(length(nPts-1,1)) ;
    if strcmpi(disp_Mode,'absolute')
        for i=1:(nPts-1)
            L0(i) = (Pts(i,:) * Pts(i+1,:)') ^ 0.5 ;
        end
    else 
        for i=1:(nPts-1)
            L0(i) = (PtsMovRef(i,:) * PtsMovRef(i+1,:)') ^ 0.5 ;
        end
    end
        
% If it is the first frame
    if frame == 1
        obj.Strains = zeros(size(Pts)) ;
    end
    
% ELSE, Compute Strains Calculation using displacement ether eulerian (relative displacement) or
% Lagrangian (absolute displacement) 

disp(['Computation of Strain ' num2str(disp_Mode)]);
if frame > 1

    % Compute Strains
    % Init
        L(1:nPts-1,1) = NaN ;
        E(1:nPts-1,1) = NaN ;
        for p = 1:(nPts-1)
            if isnan(U(p+1,1)) || isnan(U(p,1))  ; continue ; end ;
        % Actual length
            L(p) = (PtsMov(p+1,:) * PtsMov(p,:)') ^ .5 ;
        % Strains
            E(p) = (L(p)-L0(p)) / L0(p) ;
        end
        toc(t)
    % Save Strains
    obj.Strains(:,frame) = E(:) ;
end

% Display Object
    disp(obj)