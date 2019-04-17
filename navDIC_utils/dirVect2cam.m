function t3D = dirVect2cam(PtsPix , Cam) 

nbIm = length(PtsPix(1,1,:,1)) ;

% Etablissement de la direction 3D du fil a partir des deux equations des cameras
tn = zeros(nbIm,3,2) ;
for j = 1:length(Cam)
    t2 = cross(  [ squeeze(PtsPix(1,:,:,j))' - repmat(Cam(j).PPL',[nbIm 1]), repmat(Cam(j).fpix(1),[nbIm 1]) ],...
            [ squeeze(PtsPix(2,:,:,j))' - repmat(Cam(j).PPL',[nbIm 1]), repmat(Cam(j).fpix(1),[nbIm 1]) ] ) ;
    tn(:,:,j) = t2 ./ repmat(sum( t2.^2,2 ).^.5,[1 3]) ;
end
t3D = zeros(nbIm,3,2) ;

% Vecteur directeur deduit des vues cameras
for i = 1:nbIm
    t3D(i,:,3) =  cross( Cam(1).P \ tn(i,:,1)', Cam(2).P \ tn(i,:,2)' )' ;
end

for j = 1:length(Cam)
        t3D(:,:,1) =  ( Cam(1).P * t3D(:,:,3)' )' ; 
        t3D(:,:,2) =  ( Cam(2).P * t3D(:,:,3)' )' ; 
end