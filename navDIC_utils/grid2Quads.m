global hd
clearvars -except hd

Seed = hd.Seeds(2) ;

nX = 9 ;
nY = 8 ;
firstOrientation = 'x' ;

nQuads = (nX-1)*(nY-1) ;

switch firstOrientation
    case 'x'
        p1 = ( (1:nX-1) + (0:nY-2)'*nX )' ;
        p2 = p1+1 ;
        p4 = p1+nX ;
        p3 = p4+1 ;
end
Seed.Elems = [p1(:) p2(:) p3(:) p4(:)] ;
