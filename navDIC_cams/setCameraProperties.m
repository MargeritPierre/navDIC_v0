function hd = setCameraProperties(hd,num)

disp('------------------------------------------------------------------') ;
disp(['Setting of position and properties of the camera ', num2str(num),' : ',hd.Cameras(num).Name]) ;
disp('------------------------------------------------------------------') ;

if ~isempty(hd.Images)
    rep = inputdlg('utiliser image de calibration differente (oui/non)') ;
    if strcmpi(rep{1},'oui') 
        [file,path] = uigetfile('*.tif',['selectionner l''image de calibration de la camera : ',...
            hd.Cameras(num).Name,' ? utiliser. ']) ;
        img = imread([path,file]) ;
    else
        if iscell(hd.Images{1}{num})
            img = hd.Images{1}{num}{1} ;
        else
            img = hd.Images{1}{num} ;
        end
    end
else
    hd = startAllCameras(hd) ;
    img = im2single(getsnapshot(hd.Cameras(num).VidObj)) ;
end
drawToolH = drawingTool('drawROI',true ...
       ,'background', img,'title','DrawingTool : Define two points to calculate scale') ;   

if strcmpi(drawToolH.Geometries(1).Class, 'impoint')
   for p = 1:length(drawToolH.Geometries)
        pts(p,:) = obj.drawToolH.Geometries(p).Position ;
   end
elseif strcmpi(drawToolH.Geometries(1).Class, 'imline')
   pts = drawToolH.Geometries(1).Position(:,:) ;
end
refpix = sqrt(sum(diff(pts,1,1).^2,2)) ;
prompt = {['Entrer la dimension reel de la r?f?rence choisie (', num2str(refpix),' pixels) (mm) : '],...
    'Entrer la valeur de l''angle du plan horizontale de la camera 0 : Cam 1 90 : Cam 2 (?) : ',...
    'Entrer la valeur de l''angle autour de l''axe optique de la camera 0 (?) : ',...
    'Entrer le nombre de pixel dans la largeur : ',...
    'Entrer le nombre de pixel dans la hauteur : ',...
    'Entrer la taille d''un pixel: '};
title = ['Setting of position and properties of the camera ', num2str(num),' : ',hd.Cameras(num).Name];
dims = [1 100];
definput = {'10','0','0','2048','2048','0.0074'};
Camprop = inputdlg(prompt,title,dims,definput) ;

refO = str2double(Camprop{1}) ; % input('Entrer la dimension reel d''une reference du plan objet (mm) : ') ;
 % input('Entrer la dimension en pixel correspondant dans le plan image (pix) : ') ;
pixObjratio = refpix/refO ;
ang = str2double(Camprop{2}) ; % input('Entrer la valeur de l''angle du plan horizontale de la camera 0 : Cam 1 90 : Cam 2 (?) : ') ;
angTor = str2double(Camprop{3}) ; % input('Entrer la valeur de l''angle autour de l''axe optique de la camera 0 (?) : ') ;
resX = str2double(Camprop{4}) ; % input('Entrer le nombre de pixel dans la largeur : ') ;
resY = str2double(Camprop{5}) ; % input('Entrer le nombre de pixel dans la hauteur : ') ;
ix = str2double(Camprop{6})*resX ; % input('Entrer la largeur du capteur CCD : ') ;
iy = str2double(Camprop{6})*resY ; % input('Entrer la hauteur du capteur CCD : ') ;

f = 120 ; % inputdlg('Entrer la valeur de la longueur focale de l''objectif (mm) : ') ;

di = ( refpix / resY * iy + refO ) * f / refO ;
do = di * refO / ( refpix / resY * iy ) ;

angRad = ang / 180 * pi ; 
C1 = [do * sin(angRad); - do * cos(angRad); 0] ;

Ptvis1 =  [ 0; 0; 0 ] ;

angTorRad = angTor / 180 * pi ; 

Ncam = [cos(angRad) -sin(angRad) 0  ;sin(angRad) cos(angRad) 0; 0 0 1 ] * [sin(angTorRad); 0;  cos(angTorRad)] ;

Zcam = ( Ptvis1 - C1 ) / norm( Ptvis1 - C1 ) ; 
Xcam = cross( Zcam, Ncam ) ; 
Ycam = cross( Zcam, Xcam ) ; 

% sauvegarde dans la base de donnee du repere de la cam 1
hd.Cameras(num).Properties.O = C1 ;
hd.Cameras(num).Properties.X = Xcam ;
hd.Cameras(num).Properties.Y = Ycam ;
hd.Cameras(num).Properties.Z = Zcam ;


hd.Cameras(num).Properties.P = [ dot(Xcam, [1;0;0]), dot(Xcam, [0;1;0]), dot(Xcam, [0;0;1]) ;...
    dot( Ycam, [1;0;0]), dot(Ycam, [0;1;0]), dot(Ycam,[0;0;1]) ;...
    dot( Zcam, [1;0;0]), dot(Zcam, [0;1;0]), dot(Zcam,[0;0;1])] ;

hd.Cameras(num).Properties.di = di ;
hd.Cameras(num).Properties.do = do ;
hd.Cameras(num).Properties.f = f ; % mm
hd.Cameras(num).Properties.px = resX ; % nombre de pixels en x 
hd.Cameras(num).Properties.py = resY ; % nombre de pixels en y 
hd.Cameras(num).Properties.PPL = [ resX/2; resY/2] ;
hd.Cameras(num).Properties.ix = ix ; % mm largeur CCD
hd.Cameras(num).Properties.iy = iy ; % mm hauteur CCD
hd.Cameras(num).Properties.pixObjratio = pixObjratio ;
di = do * f / ( do - f ) ;

hd.Cameras(num).Properties.fpix = [di / ix * resX ; di / iy * resY ] ;

