function [t3dPC,t3dPP] = dirVect2cam(PtsPix , Cam) 

nbIm = length(PtsPix(1,1,:,1)) ;

% Etablissement de la direction 3D du fil a partir des deux equations des cameras
tnPC = zeros(nbIm,3,2) ;
tnPP = zeros(nbIm,3,2) ;
for j = 1:length(Cam)
    tPC = cross(  [ squeeze(PtsPix(1,:,:,j))' - repmat(Cam(j).PPL',[nbIm 1]), repmat(Cam(j).fpix(1),[nbIm 1]) ],...
            [ squeeze(PtsPix(2,:,:,j))' - repmat(Cam(j).PPL',[nbIm 1]), repmat(Cam(j).fpix(1),[nbIm 1]) ] ) ;
    tPP = cross( [ zeros(nbIm, 2), ones(nbIm, 1) ],...
            [ squeeze(PtsPix(2,:,:,j))' - squeeze(PtsPix(1,:,:,j))', zeros(nbIm, 1) ] ) ;
    tnPC(:,:,j) = tPC ./ repmat( sum( tPC.^2,2 ).^.5, [1 3] ) ;
    tnPP(:,:,j) = tPP ./ repmat( sum( tPP.^2,2 ).^.5, [1 3] ) ;
end
tPC = zeros(nbIm,3,2) ;
tPP = zeros(nbIm,3,2) ;

% Vecteur directeur deduit des vues cameras
for i = 1:nbIm
    t3dPC(i,:,3) =  cross( Cam(1).P \ tnPC(i,:,1)', Cam(2).P \ tnPC(i,:,2)' )' / norm(cross( Cam(1).P \ tnPC(i,:,1)', Cam(2).P \ tnPC(i,:,2)' )') ;
    t3dPP(i,:,3) =  cross( Cam(1).P \ tnPP(i,:,1)', Cam(2).P \ tnPP(i,:,2)' )' / norm(cross( Cam(1).P \ tnPP(i,:,1)', Cam(2).P \ tnPP(i,:,2)' )') ;
end

for j = 1:length(Cam)
        t3dPC(:,:,j) =  ( Cam(j).P * t3dPC(:,:,3)' )' ; 
        t3dPP(:,:,j) =  ( Cam(j).P * t3dPP(:,:,3)' )' ; 
end
end