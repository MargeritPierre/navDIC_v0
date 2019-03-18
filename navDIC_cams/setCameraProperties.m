function hd = setCameraProperties(hd,num)

f = 120 ; % inputdlg('Entrer la valeur de la longueur focale de l''objectif (mm) : ') ;
refO = input('Entrer la dimension reel d''une reference du plan objet (mm) : ') ;
refpix = input('Entrer la dimension en pixel correspondant dans le plan image (pix) : ') ;
pixel3Dratio = refpix/refO ;
ang = input('Entrer la valeur de l''angle du plan horizontalede la premiere camera 0 : Cam 1 90 : Cam 2 (°) : ') ;
resX = input('Entrer le nombre de pixel dans la largeur : ') ;
resY = input('Entrer le nombre de pixel dans la hauteur : ') ;
ix = input('Entrer la largeur du capteur CCD : ') ;
iy = input('Entrer la hauteur du capteur CCD : ') ;


di = ( refpix / resY * iy + refO ) * f / refO ;
do = di * refO / ( refpix / resY * iy ) ;

angRad = ang / 180 * pi ; 
C1 = [do * sin(angRad); - do * cos(angRad); 0] ;

Ptvis1 =  [ 0; 0; 0 ] ;

Zcam1 = ( Ptvis1 - C1 ) / norm( Ptvis1 - C1 ) ; 
Xcam1 = cross( Zcam1, [0;0;1] ) ; 
Ycam1 = cross( Zcam1, Xcam1 ) ;

% sauvegarde dans la base de donnee du repere de la cam 1
hd.Cameras(num).Properties.O = C1 ;
hd.Cameras(num).Properties.X = Xcam1 ;
hd.Cameras(num).Properties.Y = Ycam1 ;
hd.Cameras(num).Properties.Z = Zcam1 ;


hd.Cameras(num).Properties.Cam(1).P = [ dot(Xcam1, [1;0;0]), dot(Xcam1, [0;1;0]), dot(Xcam1, [1;0;0]) ;...
    dot( Ycam1, [1;0;0]), dot(Ycam1, [0;1;0]), dot(Ycam1,[1;0;0]) ;...
    dot( Zcam1, [1;0;0]), dot(Zcam1, [0;1;0]), dot(Zcam1,[1;0;0])] ;

hd.Cameras(num).Properties.do = do ;
hd.Cameras(num).Properties.f = f ; % mm
hd.Cameras(num).Properties.px = resX ; % nombre de pixels en x 
hd.Cameras(num).Properties.py = resY ; % nombre de pixels en y 
hd.Cameras(num).Properties.ix = ix ; % mm largeur CCD
hd.Cameras(num).Properties.iy = iy ; % mm hauteur CCD

di = do * f / ( do - f ) ;

hd.Cameras(num).Properties.fx = di / ix * resX ;
hd.Cameras(num).Properties.fy = di / iy * resY ;
