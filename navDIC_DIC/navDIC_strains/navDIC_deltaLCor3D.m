function [obj,hd] = navDIC_deltaLCor3D(obj,hd)

disp('deltaL_LCor3D')

nbCam = size(obj.Points,3) ; 
cocam = [2, 1] ;
sig = [dot(hd.Cameras(2).Properties.X,hd.Cameras(1).Properties.Z), dot(hd.Cameras(1).Properties.X,hd.Cameras(2).Properties.Z)] ;

refPts = camsTo3d(hd,obj.Points) ;
movPts2D(:,:,:) = obj.MovingPoints(:,:,hd.CurrentFrame,:) ;
movPts = camsTo3d(hd,movPts2D) ;

L0(1) = sqrt(sum(diff(refPts(:,:,1),1,1).^2,2)) ;
L0(2) = sqrt(sum(diff(refPts(:,:,2),1,1).^2,2)) ;
for i = 1:nbCam
    dz = sig(i) * ( movPts(:,1,cocam(i)) - refPts(:,1,cocam(i)) ) ;
    movPtsCor(:,:) = movPts(:,:,i) ;
    movPtsCor(:,:) = ( ones(size(movPtsCor(:,:))) + repmat(dz,[1 1]) / hd.Cameras(i).Properties.do ) .* movPtsCor(:,:) ;  
    L = sqrt(sum(diff(movPtsCor,1,1).^2,2)) ;
    obj.Strains(:,1,hd.CurrentFrame,i) = (L-L0(i))./L0(i) ;
end
end
% disp('deltaL_L')
%     
% % Retrieve Infos
%     frame = hd.CurrentFrame ;
%     camID = obj.CamIDs ;
%     PtsMov = round(obj.MovingPoints(:,:,frame)) ; % 
%     disp_Mode = obj.displMode ;
%     U = obj.Displacements ;
%     Pts = obj.Points; % ref absolu dl/l0
%     nPts = length(Pts(:,:)) ;
%     L0 = ones(nPts-1,1) ;
%     if strcmpi(disp_Mode,'abs') || frame==1 
%         for i=1:(nPts-1)
%             L0(i) = (Pts(i,:) * Pts(i+1,:)') ^ 0.5 ;
%         end
%     else
%         PtsMovRef = round(obj.MovingPoints(:,:,frame-1)); % ref relatif dl/l 
%         for i=1:(nPts-1)
%             L0(i) = (PtsMovRef(i,:) * PtsMovRef(i+1,:)') ^ 0.5 ;
%         end
%     end
%         
% % If it is the first frame
%     if frame == 1
%         obj.Strains = zeros(nPts-1,1) ;
%     end
%     
% % ELSE, Compute Strains Calculation using displacement ether eulerian (relative displacement) or
% % Lagrangian (absolute displacement) 
% 
% disp(['Computation of Strain ' num2str(disp_Mode)]);
% if frame > 1
% 
%     % Compute Strains
%     % Init
%         L(1:nPts-1,1) = NaN ;
%         E(1:nPts-1,1) = NaN ;
%         for p = 1:(nPts-1)
%             if isnan(U(p+1,1)) || isnan(U(p,1))  ; continue ; end ;
%         % Actual length
%             L(p) = (PtsMov(p+1,:) * PtsMov(p,:)') ^ .5 ;
%         % Strains
%             E(p) = (L(p)-L0(p)) / L0(p) ;
%         end
%     % Save Strains
%     obj.Strains(1,:,frame) = E(:) ;
% end
% 
% % Display Object
%     disp(obj)