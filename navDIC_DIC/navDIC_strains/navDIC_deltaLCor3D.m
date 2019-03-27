function [obj,hd] = navDIC_deltaLCor3D(obj,hd)
disp('deltaL_LCor3D')
Pts = obj.Points;
frame = hd.CurrentFrame ;
% Blah blah
    if frame<=obj.RefFrame
        obj.Strains = zeros([size(Pts,1)/2 2 frame]) ;
    end
nbCam = size(obj.Points,3) ; 

cocam = [2, 1] ;

refPts = camsTo3d(hd,obj.Points) ;
movPts2D(:,:,:) = obj.MovingPoints(:,:,hd.CurrentFrame,:) ;
movPts = camsTo3d(hd,movPts2D) ;

for i = 1:nbCam
    L02D = sqrt(sum( ( obj.Points(2:2:end,:,i) - obj.Points(1:2:end,:,i) ).^2,2) ) ;
    L2D = sqrt(sum( ( movPts2D(2:2:end,:,i) - movPts2D(1:2:end,:,i) ).^2,2) ) ;
    dz = ( movPts(:,:,cocam(i)) - refPts(:,:,cocam(i)) ) *...
        [dot(hd.Cameras(cocam(i)).Properties.X,hd.Cameras(i).Properties.Z);...
        dot(hd.Cameras(cocam(i)).Properties.Y,hd.Cameras(i).Properties.Z); 0];
    movPtsCor(:,:) = ( ones(size(movPts(:,:,i))) + repmat(dz,[1 1]) / hd.Cameras(i).Properties.do ) .* movPts(:,:,i) ;  
    L0 = sqrt(sum(( refPts(2:2:end,:,i) - refPts(1:2:end,:,i) ).^2,2)) ;
    L = sqrt(sum( ( movPtsCor(2:2:end,:) - movPtsCor(1:2:end,:) ).^2,2) ) ;
    obj.Strains(:,2,hd.CurrentFrame,i) = (L2D-L02D)./L02D ; % 3D
    obj.Strains(:,1,hd.CurrentFrame,i) = (L-L0)./L0 ; % 2D
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