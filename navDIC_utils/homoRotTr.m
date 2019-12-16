function [ptsPixCor,ptsRec] = homoRotTr(ptsPix,ref,Cam,a,uz)

% ---------------------------------------------------------------
% Inputs
% ---------------------------------------------------------------
% a : liste angle hors plan (en radian)
% uz : deplacement hors plan (en millimetre)
% ref : id du point de reference
% ptsPix : coordonnées des points apparent en pixels

% ---------------------------------------------------------------
% Outputs
% ---------------------------------------------------------------
% ptsPixCor : coordonnées des points apparent en pixels
% ptsRec : coordonnées des points apparent en pixels

xc = Cam.PPL(1) ;
yc = Cam.PPL(2) ;
x1 = mean(ptsPix(ref,1),1) ;
y1 = mean(ptsPix(ref,2),1) ;
f = Cam.fpix(2) ;
nPts = size(ptsPix,1) ;

ptsPixCor = zeros( size( ptsPix ) ) ;
ptsRec = zeros( size( ptsPix, 1 ), 3 ) ;

if abs(Cam.X'*[0;0;1]) > abs(Cam.Y'*[0;0;1])
    % Applique l'homographie pour corriger la rotation hors-plan
    
    % Abscisses
    ptsPixCor(:,1) = ( ptsPix(:,1) - repmat(x1,[nPts 1]) ) * f ./ ( repmat( f * cos(a) ,[nPts 1]) - ( ptsPix(:,1) - repmat(xc,[nPts 1]) ) * sin(a) ) + repmat(x1,[nPts 1])  ;
    ptsRec(:,1) = ( ptsPixCor(:,1) - repmat(x1,[nPts 1]) ) * cos(a) + repmat(x1,[nPts 1]) ;
    
    % positions hors plan
    ptsRec(:,3) = repmat(f, [nPts 1]) + ( ptsPixCor(:,1) - repmat(x1,[nPts 1]) ) * sin(a) ;
    
    % Ordonnées
    ptsPixCor(:,2) = ( ptsPix(:,2) - repmat(yc,[nPts 1]) ) .* ( repmat( f ,[nPts 1]) + ( ptsPixCor(:,1) - repmat(x1,[nPts 1]) ) * sin(a)  ) / f + repmat(yc,[nPts 1]) ;
    ptsRec(:,2) = ptsPixCor(:,2) ;
    
else
    % Applique l'homographie pour corriger la rotation hors-plan
    % Ordonnées
    ptsPixCor(:,2) = ( ptsPix(:,2) - repmat(y1,[nPts 1]) ) * f ./ ( repmat( f * cos(a) ,[nPts 1]) - ( ptsPix(:,2) - repmat(yc,[nPts 1]) ) * sin(a) ) + repmat(y1,[nPts 1])  ;
    ptsRec(:,2) = ( ptsPixCor(:,2) - repmat(y1,[nPts 1]) ) * cos(a) + repmat(y1,[nPts 1]) ;
    
    % positions hors plan
    ptsRec(:,3) = repmat(f, [nPts 1]) + ( ptsPixCor(:,2) - repmat(y1,[nPts 1]) ) * sin(a) ;
    
    % Abscisses
    ptsPixCor(:,1) = ( ptsPix(:,1) - repmat(xc,[nPts 1]) ) .* ( repmat( f ,[nPts 1]) + ( ptsPixCor(:,2) - repmat(y1,[nPts 1]) ) * sin(a)  ) / f + repmat(xc,[nPts 1]) ;
    ptsRec(:,1) = ptsPixCor(:,1) ;
end

% Applique le scaling pour corriger la translation hors-plan
ptsPixCor = ( ptsPixCor - repmat( [xc, yc], [nPts 1]) ) *( f + uz * Cam.pixObjratio ) / f + repmat([xc, yc],[nPts 1]) ;
ptsRec = ( ptsRec - repmat( [xc, yc, 0], [nPts 1]) ) * ( f + uz * Cam.pixObjratio ) / f ;

% passage du plan image au plan objet pour les points reconstruits
ptsRec = ptsRec / Cam.pixObjratio ;
ptsRec(:,:,2) = ( Cam.P \ ptsRec(:,:)' )' + repmat( Cam.O', [nPts 1] ) ;
end



