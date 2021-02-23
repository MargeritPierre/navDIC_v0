%% Preview of a camera with custom functionalities

global hd
%% INITIALIZE THE FIGURES
img = getsnapshot(hd.Cameras.VidObj);
iii = image(img);
colormap('gray');
preview(hd.Cameras.VidObj,iii);
axis equal

%% Lecture des propriétés de la source :
src = hd.Cameras.VidObj.Source;

set(src)
%% Modification du ExposureTime:

get(src,'ExposureTimeAbs')
set(src,'ExposureTimeAbs',5000);

%% Réglage du FrameRate :
disp(['AcquisitionFrameRateAbs before : ', num2str(get(src,'AcquisitionFrameRateAbs'))]);
set(src,'AcquisitionFrameRateAbs',3);

%% [OLD] INITIALIZE THE FIGURES
%% INITIALIZE THE FIGURES
    seed = hd.Seeds(1) ;
    contourTag = 'J-integral Contours' ;
    valuesTag = 'J-integral Values' ;
    % Init the contour figure
        figContours = findobj(0,'Name',contourTag) ;
        if isempty(figContours)
            figContours = figure('Name',contourTag) ; 
            figContours.Position = [figContours.Position(1:2)+figContours.Position(3:4)/2.*[0 .5] figContours.Position(3:4)/2] ;
        end
        clf(figContours) ;
        refImg = hd.Images{1}{1} ;
        im = imagesc(repmat(refImg(:,:,1),[1 1 3])) ;
        mesh = patch('Vertices',seed.MovingPoints(:,:,1),'Faces',seed.Elems,'Facecolor','w','edgecolor','none','Facealpha',0.5) ;
        axis tight
        axis equal
        axis ij
        axis off
        set(gca,'xtick',[],'ytick',[]) ; 
        [nI,nJ] = size(refImg) ;
    % Add a ruler line
        ruler = images.roi.Line ;
        ruler.Parent = gca ;
        ruler.Color = 'r' ;
        ruler.Position = [1 1 ; nJ 1] ; %[1/20 0.095 ; 1/20 0.865].*[nJ nI] ;
        rulerLength = @(ruler)sqrt(sum(diff(ruler.Position,1,1).^2,2)) ;
        addlistener(ruler,'MovingROI',@(src,evt)disp(['length: ',num2str(rulerLength(src))])) ;
        
%% 
