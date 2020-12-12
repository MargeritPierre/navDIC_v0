%% APPLY MEDIAN FILTERING ON A SERIES OF IMAGES

%% FOR RELATIVELY SMALL IMAGE SERIES: vectorization
clearvars -except hd
global hd

cam = 1 ;
medFiltSize = 10 ;

disp('')
disp('------ Median Filtering -------')

% Retrieving the images
disp('   Concatenate images..')
IMG = cat(4,hd.Images{cam}{:}) ;

% Convert to a 2D array
disp('   to 2D array..')
[nI,nJ,nC,nIMG] = size(IMG) ;
IMG = reshape(IMG,[],nIMG) ;

% Median Filtering
disp('   Median Filtering..')
IMG = medfilt2(IMG,[1 medFiltSize],'symmetric') ;

% Reshaping back
disp('   Reshaping..')
IMG = reshape(IMG,[nI nJ nC nIMG]) ;
IMG = num2cell(IMG,1:3) ;


%% FOR LONG IMAGE SERIES (Huge data): Loop (not implemented yet)
clearvars -except hd
global hd

cam = 1 ;
medFiltSize = 5 ;

disp('')
disp('------ Median Filtering -------')

% Retrieving the images
disp('   Retrieve images..')
IMG = hd.Images{cam} ;

%% PUSH TO navDIC
newCam = 3 ;
namePrefix = ['Median ' num2str(medFiltSize)] ;

disp('   Pushing to navDIC..')
hd.Cameras(newCam) = hd.Cameras(cam) ;
hd.Cameras(newCam).Name = [namePrefix ' | ' hd.Cameras(cam).Name] ;
hd.Images{newCam} = IMG(:) ;

