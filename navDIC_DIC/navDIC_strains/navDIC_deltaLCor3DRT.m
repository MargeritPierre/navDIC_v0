function [obj,hd] = navDIC_deltaLCor3DRT(obj,hd)
disp('deltaL_LCor3DRT')
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
    Proj = [dot(hd.Cameras(cocam(i)).Properties.X,hd.Cameras(i).Properties.Z) ; ...
        dot(hd.Cameras(cocam(i)).Properties.Y,hd.Cameras(i).Properties.Z); 0] ;
    dzR0 = meanNoNaN( ( refPts(2:2:end,:,cocam(i)) - refPts(1:2:end,:,cocam(i)) ) * Proj  , 1) ;

    L02D = sqrt(sum( ( obj.Points(2:2:end,:,i) - obj.Points(1:2:end,:,i) ).^2,2) ) ;
    L2D = sqrt(sum( ( movPts2D(2:2:end,:,i) - movPts2D(1:2:end,:,i) ).^2,2) ) ;
    dz = ( movPts(:,:,cocam(i)) - refPts(:,:,cocam(i)) ) * Proj ;
    
    movPtsCor(:,:) = ( ones(size(movPts(:,:,i))) + ( repmat(dz,[1 size(movPts(:,:,i),2)]) +...
        repmat([0; dzR0],[size(movPts(:,:,i),1)/2 size(movPts(:,:,i),2)]) ) /...
        hd.Cameras(i).Properties.do ) .* movPts(:,:,i) ; 
    refPtsCor(:,:) = ( ones(size(refPts(:,:,i))) +...
        repmat([0; dzR0],[size(refPts(:,:,i),1)/2 size(refPts(:,:,i),2)]) /...
        hd.Cameras(i).Properties.do ) .* refPts(:,:,i) ; 
    L0 = sqrt(sum(( refPtsCor(2:2:end,:) - refPtsCor(1:2:end,:) ).^2,2)) ;
    L = sqrt(sum( ( movPtsCor(2:2:end,:) - movPtsCor(1:2:end,:) ).^2,2) ) ;
    obj.Strains(:,2,hd.CurrentFrame,i) = (L2D-L02D)./L02D ; % 3D
    obj.Strains(:,1,hd.CurrentFrame,i) = (L-L0)./L0 ; % 2D
end

end