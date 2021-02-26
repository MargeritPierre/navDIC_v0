%% Preview of a camera with custom functionalities

global hd
%% INITIALIZE THE FIGURES
fig= figure('Name','Preview custom'); 
fig.Position = [fig.Position(1:2)+fig.Position(3:4)/2.*[0 .5] fig.Position(3:4)/2] ;

refImg =getsnapshot(hd.Cameras.VidObj) ;
iii = image(refImg);
colormap('gray');
axis equal;
%% Preview à lancer dans un second temps. Pourquoi pas dans le même bloc ? raison inconnue!
preview(hd.Cameras.VidObj,iii);

%% Lecture des propriétés de la source :
src = hd.Cameras.VidObj.Source;
set(src)

%% Modification du ExposureTime:
disp(['ExposureTimeAbs before : ', num2str(get(src,'ExposureTimeAbs'))]);

set(src,'ExposureTimeAbs',200000);

%% Réglage du FrameRate :
disp(['AcquisitionFrameRateAbs before : ', num2str(get(src,'AcquisitionFrameRateAbs'))]);
set(src,'AcquisitionFrameRateAbs',3);


%% Ajout d'une ligne 
[nI,nJ] = size(refImg);
% Add a ruler line
ruler = images.roi.Line ;
ruler.Parent = gca ;
ruler.Color = 'r' ;
ruler.Label = 'ruler 1';
ruler.MarkerSize=4;
ruler.Position = [1 1 ; nJ 1] ; %[1/20 0.095 ; 1/20 0.865].*[nJ nI] ;
rulerLength = @(ruler)sqrt(sum(diff(ruler.Position,1,1).^2,2)) ;
rulerAngle = @(ruler)atan2( ruler.Position(2,2) -ruler.Position(1,2) , ruler.Position(2,1) -ruler.Position(1,1)) * 180 / pi;

addlistener(ruler,'MovingROI',@(src,evt)disp([ruler.Label, ' –-  length: ',num2str(rulerLength(src))])) ;
addlistener(ruler,'MovingROI',@(src,evt)disp([ruler.Label, ' --  angle: ',num2str(rulerAngle(src))])) ;


%% Ajout de deux lignes avec symétrie par rapport au plan médian vertical
[nI,nJ] = size(refImg);
% Add a ruler line
ruler_master = images.roi.Line('Parent', gca, 'Color','r','Label','master','MarkerSize',4) ;
ruler_slave = images.roi.Line('Parent', gca, 'Color','g','Label','slave','MarkerSize',4) ;

ruler_master.Position = [1 1 ; nJ/2 1] ;
ruler_slave.Position = [nJ/2 + (nJ/2 - ruler_master.Position(1,1)) ruler_master.Position(1,2);  nJ/2 + (nJ/2 - ruler_master.Position(2,1)) ruler_master.Position(2,2)];

adjust_pos_slave = @(master)set(ruler_slave,'Position',[nJ/2 + (nJ/2 - master.Position(1,1)) master.Position(1,2);  nJ/2 + (nJ/2 - master.Position(2,1)) master.Position(2,2)]);


addlistener(ruler_master,'MovingROI',@(src,evt)adjust_pos_slave(src));


%% Ajout d'une croix avec les axes horizontal/vertical
[nI,nJ] = size(refImg);
% Add a ruler line
crosshair = images.roi.Crosshair('Parent', gca, 'Color','b') ;
crosshair.Position =  [nJ/2,nI/2];