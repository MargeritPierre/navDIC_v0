%% OCR ON NAVDIC IMAGES
global hd

seed = hd.Seeds(end) ; % Seed for ROI
frames = 1:hd.nFrames ; hd.CurrentFrame ;

bboxRatio = 1.*[.5 .8] ;
subSampling = 1 ;
blackandwhite = true ;
normalize = true ;
backgroundFiltSz = 42.*[1 1] ;

% Text bounding boxes (given by the seed position)
    bbox = seed.MovingPoints(:,:,frames) ;
    v = (bbox(1,:,:)-bbox([2 4],:,:)).*(1-bboxRatio(:)) ;
    bbox = bbox + v(1,:,:).*[-1;1;1;-1]/2 + v(2,:,:).*[-1;-1;1;1]/2 ;
% Pixel coordinates for interpolation
    Ld = round(sqrt(sum((bbox(1,:,:)-bbox([2 4],:,:)).^2,2)))*subSampling ;
    Lmax = max(Ld,3) ;
    [e1,e2] = meshgrid((0:Lmax(1))/Lmax(1),(1:Lmax(2))/Lmax(2)) ;
    x = bbox(1,:,:).*(1-e1(:)).*(1-e2(:)) ...
            + bbox(2,:,:).*(e1(:)).*(1-e2(:)) ...
            + bbox(3,:,:).*(e1(:)).*(e2(:)) ...
            + bbox(4,:,:).*(1-e1(:)).*(e2(:)) ;
% Image processing
    IMG = cell(size(hd.Images{seed.CamIDs})) ;
    wtbr = waitbar(0,'Image processing..') ;
    for ff = 1:numel(frames)
    % Get original
        fr = frames(ff) ;
        img0 = hd.Images{seed.CamIDs}{fr} ;
    % Black and white
        if blackandwhite ; img0 = sum(img0/size(img0,3),3,'native') ; end
    % Interpolation
        IMG{fr} = interp2(single(img0(:,:,1)),x(:,1,ff),x(:,2,ff),'linear',0) ;
        for cc = 2:size(img0,3) 
            IMG{fr}(:,cc) = interp2(single(img0(:,:,cc)),x(:,1,ff),x(:,2,ff),'linear',0) ;
        end
    % Reshape
        IMG{fr} = reshape(IMG{fr},size(e1,1),size(e1,2),[]) ;
    % Filtering for background removal
        IMG{fr} = abs(IMG{fr}-medfilt2(IMG{fr},backgroundFiltSz,'symmetric')) ;
    % Recast
        %IMG{fr} = cast(IMG{fr},class(img0)) ;
        IMG{fr} = flip(IMG{fr},1) ;
    % Normalize
        if normalize ; IMG{fr} = IMG{fr}*(max(getrangefromclass(IMG{fr}))/max(IMG{fr}(:))) ; end
    % Wait bar
        wtbr = waitbar(ff/numel(frames),wtbr) ;
    end
    delete(wtbr) ;
    
%% PUSH TO NAVDIC
    newCamIdx = 3 ; numel(hd.Cameras)+1 ;
    newCam = hd.Cameras(seed.CamIDs) ;
    newCam.Name = ['Process Text | ' hd.Cameras(seed.CamIDs).Name] ;
    hd.Images{newCamIdx} = IMG ;
    hd.Cameras(newCamIdx) = newCam ;


%% BUILTIN MATLAB OCR
frames = hd.CurrentFrame ;
for fr = frames(:)'
% Prepare image
    img = IMG{fr} ;
% Text recognition
    ocrResults = ocr(img,'TextLayout','Word','CharacterSet','0123456789') 
    if ~isempty(ocrResults.Text)
        img = insertObjectAnnotation(img,'rectangle',ocrResults.WordBoundingBoxes,ocrResults.WordConfidences);
    elseif size(img,3)==1
        img = repmat(img,1,1,3) ;
    end
    %figure ; 
    clf ; image(img) ; axis equal ij
end

%% CUSTOM DETERMINISTIC OCR FOR DIGITAL NUMBERS

bbox = [0 4 ; 80 33] ;
nDigits = 5 ;
digitPow = 10.^(3:-1:-1) ;
frames = 1:hd.nFrames ; hd.CurrentFrame+4 ;
maskShape = 'quad' ; % 'tri' or 'quad' ;
imgTrsh = 45/100 ; imgAmpl = 12 ; % sigmoid thresholding
maskTrsh = .25 ;

% Define masks for each digit
Ld = diff(bbox,1)./[nDigits 1] ;
digitBoxes = bbox(1,:)+[0;1]*Ld + reshape([0:nDigits-1],1,1,nDigits).*[1 0].*Ld ;
nodes = [digitBoxes([1 2 2 1],1,:) digitBoxes([1 1 2 2],2,:)] ;
nodes = [nodes ; mean(nodes([2 3],:,:),1) ; mean(nodes([1 4],:,:),1)] ;

% Image coordinates
img0 = IMG{frames(1)} ;
[xx,yy] = meshgrid(1:size(img0,2),1:size(img0,1)) ;

switch maskShape
    case 'tri'
    % Define triangles for each mask
        nodes = [nodes ; mean(nodes([1 2 5 6],:,:),1) ; mean(nodes([3 4 5 6],:,:),1)] ;
        tri = [1 2 7 ; 2 5 7 ; 5 6 7 ; 6 1 7 ; 5 8 6 ; 5 3 8 ; 3 4 8 ; 4 6 8] ;
        poly = reshape(nodes(tri,:,:),8,3,2,nDigits) ;
        poly(:,end+1,:,:) = NaN ;
        poly(end+1,:,:,:) = nodes([6 7 5 8],:,:) ;
    % Digit correspondence
        digitCondition = logical([ ...
            0 0 0 0 0 0 0 0 0 ; ... .
            1 1 0 1 0 1 1 1 0 ; ... 0
            0 1 0 0 0 1 0 0 0 ; ... 1
            1 1 0 0 0 0 1 1 1 ; ... 2
            1 1 0 0 0 1 1 0 1 ; ... 3
            0 1 0 1 0 1 0 0 1 ; ... 4
            1 0 0 1 0 1 1 0 1 ; ... 5
            1 0 0 1 0 1 1 1 1 ; ... 6
            1 1 0 0 0 1 0 0 0 ; ... 7
            1 1 0 1 0 1 1 1 1 ; ... 8
            1 1 0 1 0 1 1 0 1 ; ... 9
            ]) ;
    case 'quad'
        side = 6 ;
    % Define triangles for each mask
        nodes = reshape(nodes([1 2 , 2 5 , 5 3 , 3 4 , 4 6 , 6 1 , 5 6],:,:),2,7,2,nDigits) ;
        nodes = reshape(mean(nodes,1),7,2,nDigits) ;
        nodes = nodes + [0 1;-1 0;-1 0;0 -1;1 0;1 0; 0 0]*side/2 ;
        poly = reshape(nodes,7,1,2,nDigits) + reshape([-1 -1;1 -1;1 1;-1 1],1,4,2)*side/2 ;
    % Digit correspondence
        digitCondition = logical([ ...
            0 0 0 0 0 0 0 ; ... .
            1 1 1 1 1 1 0 ; ... 0
            0 1 1 0 0 0 0 ; ... 1
            1 1 0 1 1 0 1 ; ... 2
            1 1 1 1 0 0 1 ; ... 3
            0 1 1 0 0 1 1 ; ... 4
            1 0 1 1 0 1 1 ; ... 5
            1 0 1 1 1 1 1 ; ... 6
            1 1 1 0 0 0 0 ; ... 7
            1 1 1 1 1 1 1 ; ... 8
            1 1 1 1 0 1 1 ; ... 9
            ]) ;
end

% Create masks
nMasksByDigit = size(poly,1) ;
nPolyNodes = size(poly,2) ;
X = reshape(permute(poly,[4 1 2 3]),nDigits*nMasksByDigit,nPolyNodes,2) ;
masks = NaN(numel(img0),nDigits,nMasksByDigit) ;
for pp = 1:size(X,1)
    masks(inpolygon(xx(:),yy(:),X(pp,:,1),X(pp,:,2)),pp) = 1 ;
end
nCond = max(eps,sum(digitCondition,2)) ;


% Init figure
clf ; im = imagesc(repmat(img0,[1 1 3])) ; axis equal ij ; %colormap gray
pa = patch(X(:,:,1)',X(:,:,2)',X(:,:,1)'*0,'FaceColor','flat','EdgeColor','r','FaceAlpha',.7) ;
pa.FaceVertexCData = NaN(nMasksByDigit*nDigits,1) ;

% Detect numbers
numbers = NaN(hd.nFrames,1) ;
for fr = frames(:)'
    img = IMG{fr} ;
    img = 1./(1+exp(-imgAmpl*(img-imgTrsh))) ;
    img = (img-min(img(:)))/range(img(:)) ;
    data = img(:).*masks ; % [nPixels nDigits nMasksByDigit]
    mData = mean(data,1,'omitnan') ; % [1 nDigits nMasksByDigit]
    pa.FaceVertexCData = reshape(mData(:),[],1) ;
    %pa.FaceVertexCData(pa.FaceVertexCData<.2) = NaN ;
    mData = permute(mData,[1 3 2]) ; % [1 nMasksByDigit nDigits]
    confidence = sum((mData-maskTrsh).*(digitCondition-.5),2) ;
%     confidence = sum(digitCondition.*mData,2)./nCond ...
%                     - sum((~digitCondition).*mData,2)./(9-nCond) ;
    [c,n] = max(confidence(:,:),[],1) ;
    num = max((n-2),0) , c
    if ~any(isnan(c)) ; numbers(fr) = sum(num.*digitPow) ; end
    im.CData = repmat(img,[1 1 3]) ; drawnow ;
end



%%

imgTrsh = 50/100 ;
imgAmpl = 20 ; 
x = linspace(0,1,1000) ;
f = 1./(1+exp(-imgAmpl*(x-imgTrsh))) ;

clf ; plot(x,f) ;





