function PtsPixCor = PixCor2cam(PtsPix, Cam, imRef) 
nbIm = size(PtsPix,3);
cocam = [2 1] ;
PtsPixCor = zeros( size(PtsPix,1), 3, nbIm, length(Cam) ) ;
t3D = dirVect2cam(PtsPix , Cam) ;

for i = 1:length(Cam)
% Camera plus correction deplacement 
    % 1- calcul mouvement hors plan
    dzRef(:,1, i) = ( squeeze(PtsPix(1,:,:,cocam(i)))' - repmat(PtsPix(1,:,imRef,cocam(i)),[nbIm 1]) ) /...
        Cam(cocam(i)).fpix(1) * Cam(cocam(i)).do *...
        round([dot(Cam(cocam(i)).X,Cam(i).Z); dot(Cam(cocam(i)).Y,Cam(i).Z)]) / Cam(i).do * Cam(i).fpix(1) ;
    % 2- recalcul position points corrigees
    R(:,:,i) = [ squeeze(PtsPix(1,:,:,i))' - repmat(Cam(i).PPL',[nbIm 1]),repmat(Cam(i).fpix(1),[nbIm 1])] .*...
        repmat( ( 1 + dzRef(:,1, i) / Cam(i).fpix(1) ),[1,3] ) ;
    for j = 1:nbIm
        PtsPixCor(1,:,j,i) = R(j,:,i) ;
        X = [ [ (PtsPix(2,:,j,i)-Cam(i).PPL'), Cam(i).fpix(1)]' , -t3D(j,:,i)' ] ;
        A = X \ R(j,:,i)' ;
        PtsPixCor(2,:,j,i) = A(1) * [(PtsPix(2,:,j,i)-Cam(i).PPL'), Cam(i).fpix(1)] ;
    end

end


end