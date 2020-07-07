%% APPLY MEDIAN FILTERING ON A SERIES OF IMAGES
clearvars -except hd
global hd

cam = 1 ;
medFiltSize = 3 ;
newCam = 3 ;
newCamName = 'MedianFiltered' ;

disp('')
disp('------ Median Filtering -------')

% Retrieving the images
disp('   Retrieve images..')
IMG = hd.Images{cam} ;
IMG = cat(4,IMG{:}) ;

% Convert to a 2D array
disp('   to 2D array..')
sz = size(IMG) ;
IMG = reshape(IMG,[],sz(end)) ;

% Median Filtering
img = IMG ;
disp('   Median Filtering..')
wtbr = waitbar(0,'Median Filtering...') ;
for fr = 1:size(img,2)
    ii = fr-floor(medFiltSize/2) + (0:medFiltSize-1) ;
    ii = min(ii,size(img,2)) ;
    ii = max(ii,1) ;
    img(:,fr) = median(img(:,ii),2) ;
    wtbr = waitbar(fr/size(img,2),wtbr,['Median Filtering... (' num2str(fr) '/' num2str(size(img,2)) ')']) ;
    drawnow ;
end
delete(wtbr) ;

% Reshaping back
disp('   Reshaping..')
img = reshape(img,sz) ;
img = num2cell(img,1:3) ;

% Pushing to navDIC
disp('   Pushing to navDIC..')
hd.Cameras(newCam) = hd.Cameras(cam) ;
hd.Cameras(newCam).Name = newCamName ;
hd.Images{newCam} = img(:) ;

