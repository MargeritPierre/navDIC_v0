%% INITALIZE GLOBAL DIC

% TIME THE CODE ?
    profile off
    if codeProfile; profile on ; end

% RETRIEVE THE DATA
    navDICFrames = 1:hd.nFrames ;
    if ischar(frames) ; frames = eval(['navDICFrames(',frames,')']) ; end
    IMG = cat(4,hd.Images{camID}{frames}) ; % /!\ STILL NO SUPPORT FOR MULTICOLOR IMAGES YET
 
% PROCESS IMAGES
    [nI,nJ,nFrames] = size(IMG) ;
    [JJ,II] = meshgrid(1:nJ,1:nI) ;
    % Img Normalization
        imgClassRange = 1 ;
        if normToImageClassRange
            imgClassRange = double(range(getrangefromclass(IMG(1)))) ;
        end
    
% TIME MEANING
    if timeMeanLength>0
        meanTime = 2*timeMeanLength+1 ;
        TimeKernel = ones(meanTime,1)/meanTime ;
        IMG0 = IMG ;
        IMG = IMG*0 ;
        wtbr = waitbar(0,'Time Averaging') ;
        for ii = 1:nFrames
            for tt = -timeMeanLength:timeMeanLength
                ind = max(1,min(nFrames,ii+tt)) ;
                IMG(:,:,:,ii) = IMG(:,:,:,ii) + IMG0(:,:,:,ind)*TimeKernel(tt+timeMeanLength+1) ;
            end
            wtbr = waitbar(ii/nFrames,wtbr,['Time Averaging (',num2str(ii),'/',num2str(nFrames),')']) ;
        end
        delete(wtbr) ;
        clear IMG0
    end

% REFERENCE IMAGE
    % Correct reference frale index
        switch refFrame
            case 'first'
                refFrame = 1 ;
            case 'last'
                refFrame = size(IMG,4) ;
            otherwise % A number has been given
        end
    % Indices of non-used/average frames
        switch dicDir
            case 1
                avgFrames = 1:refFrame ;
                dicFrames = refFrame+1:nFrames ;
            case -1
                avgFrames = refFrame:nFrames ;
                dicFrames = refFrame-1:-1:1 ;
        end
    img0 = IMG(:,:,:,refFrame) ;
    if averagePreviousFrames % AVERAGING THE FIRST OR LAST IMAGES
        img0 = IMG(:,:,:,1)*0 ;
        for ii = avgFrames
            img0 = img0 + IMG(:,:,:,ii)/length(avgFrames) ;
        end
    end
    
% INIT THE FIGURE
    figDIC = findobj(groot,'tag',figTag) ;
        if isempty(figDIC)
            figDIC = figure('tag',figTag,'name',figTag,'numbertitle','off') ;
        end
    
    