%% SOME EXAMPLE LINES FOR USING navDIC
global hd
time = hd.TimeLine ;

%% ACESS IMAGES
cameraID = 1 ;
frame = hd.CurrentFrame ;
nFrames = hd.nFrames ;

img = hd.Images{cameraID}{frame} ;

figure ; 
axis equal tight ij ; 
imagesc(img)
colormap gray

%% ACCESS DIC RESULTS

seed = hd.Seeds(1);

refImg = seed.refImgs ; % reference image(s)
X = seed.Points ; % reference config [nPoints, nCoord==2]
x = seed.MovingPoints ; % all configs [nPoints, nCoord==2, nFrames]
u = x-X ; % all displacements [nPoints, nCoord==2, nFrames]
time = sum(hd.TimeLine.*cumprod([0 0 0 24 60 60],'reverse'),2) ; % [nFrames 1]

figure ; 
xplot = permute(x,[3 1 2]) ;
plot3(xplot(:,:,1),xplot(:,:,2),time-time(1),'o-')
xlabel 'X (pix)'
ylabel 'Y (pix)'
ylabel 'time (secs)'


%% POST-COMPUTED DATAFIELDS

seed.computeDataFields() ; % datafields may not have been comuted already..
normU = seed.DataFields.U ; % [nPoints, 1 , nFrames]

figure ; plot(time-time(1),normU(:,:).')