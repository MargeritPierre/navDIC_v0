function [obj,hd] = navDIC_deltaLCor3DRT(obj,hd)
disp('deltaL_LCor3DRT')
Pts = obj.Points;
frame = hd.CurrentFrame ;

% Blah blah
    if frame<=obj.RefFrame
        obj.Strains = zeros([size(Pts,1)/2 2 frame]) ;
    end
    
nbCam = size(obj.Points,3) ;
for i = 1:nbCam
    Cam(i) = hd.Cameras(i).Properties ;
end

   
cocam = [2, 1] ;

movPts2D = obj.MovingPoints(:,:,hd.CurrentFrame,:) ;

% Methode de correction par determination direction fil

PtsPix(:,:,1,:) = obj.Points ; % ref points
PtsPix(:,:,2,:) = obj.MovingPoints(:,:,hd.CurrentFrame,:) ; % moving points

PtsPixCor = PixCor2cam(PtsPix, Cam) ;

for i = 1:length(Cam)
    % pixels apparents
    L02D = sqrt(sum( ( obj.Points(2:2:end,:,i) - obj.Points(1:2:end,:,i) ).^2,2) ) ;
    L2D = sqrt(sum( ( obj.MovingPoints(2:2:end,:,hd.CurrentFrame,i) -...
        obj.MovingPoints(1:2:end,:,hd.CurrentFrame,i) ).^2,2) ) ;
    % pixels corrig?s
    L0 = sum( ( PtsPixCor(2:2:end,:,1,i) - PtsPixCor(1:2:end,:,1,i) ).^2, 2 ).^0.5 ;
    L = sum( ( PtsPixCor(2:2:end,:,2,i) - PtsPixCor(1:2:end,:,2,i) ).^2 ).^0.5 ;
    % deformations correspondantes
    obj.Strains(:,2,hd.CurrentFrame,i) = (L2D-L02D)./L02D ; % 2D
    obj.Strains(:,1,hd.CurrentFrame,i) = (L-L0)./L0 ; % 3D
end


% Methode de correction par points
% for i = 1:nbCam
%     Proj = [dot(hd.Cameras(cocam(i)).Properties.X,hd.Cameras(i).Properties.Z) ; ...
%         dot(hd.Cameras(cocam(i)).Properties.Y,hd.Cameras(i).Properties.Z); 0] ;
%     dzR0 = meanNoNaN( ( refPts(2:2:end,:,cocam(i)) - refPts(1:2:end,:,cocam(i)) ) * Proj  , 1) ;
% 
%     
%     dz = ( movPts(:,:,cocam(i)) - refPts(:,:,cocam(i)) ) * Proj ;
%     
%     movPtsCor(:,:) = ( ones(size(movPts(:,:,i))) + ( repmat(dz,[1 size(movPts(:,:,i),2)]) +...
%         repmat([0; dzR0],[size(movPts(:,:,i),1)/2 size(movPts(:,:,i),2)]) ) /...
%         hd.Cameras(i).Properties.do ) .* movPts(:,:,i) ; 
%     refPtsCor(:,:) = ( ones(size(refPts(:,:,i))) +...
%         repmat([0; dzR0],[size(refPts(:,:,i),1)/2 size(refPts(:,:,i),2)]) /...
%         hd.Cameras(i).Properties.do ) .* refPts(:,:,i) ; 
%     L0 = sqrt(sum(( refPtsCor(2:2:end,:) - refPtsCor(1:2:end,:) ).^2,2)) ;
%     L = sqrt(sum( ( movPtsCor(2:2:end,:) - movPtsCor(1:2:end,:) ).^2,2) ) ;
%     obj.Strains(:,2,hd.CurrentFrame,i) = (L2D-L02D)./L02D ; % 3D
%     obj.Strains(:,1,hd.CurrentFrame,i) = (L-L0)./L0 ; % 2D
% end

end