if 1 % USE THIS TO GO DIRECTLY TO DIC

clc
global hd
clearvars -except hd

% PARAMETERS
    camID = 1 ;
    seedNumber = 2 ;
    frames = '1:1:295' ; % Frames taken for DIC (allows decimation)
    dicDir = -1 ; % DIC running direction ('forward=1' or 'backward=-1')
    refFrame = 'last' ; % Reference image ('first' , 'last' or number)
    averagePreviousFrames = false ; % Ref frame is the average of the previous/next ones in forward/backward modes
    normToImageClassRange = true ; % Normalize images to their dataclass range
    timeMeanLength = 0 ; % Time averaging of images
    codeProfile = false ; % Code timing

    
% TIME THE CODE ?
    profile off
    if codeProfile; profile on ; end

% RETRIEVE THE DATA
    navDICFrames = 1:hd.nFrames ;
    frames = eval(['navDICFrames(',frames,')']) ;
    IMG = hd.Images{camID}(:,:,1,frames) ; % /!\ STILL NO SUPPORT FOR MULTICOLOR IMAGES YET
 
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
    figGlobalDIC = findobj(groot,'tag','globalDICfigure') ;
        if isempty(figGlobalDIC)
            figGlobalDIC = figure ;
        end
    
%% PLAY VIDEO OF PROCESSED IMAGES
    clf(figGlobalDIC,'reset') ;
    figure(figGlobalDIC) ;
    axes('position',[0 0 1 1])
        im = imagesc(1:nJ,1:nI,img0) ; colormap(gray)
        ttl = title('','interpreter','none','units','normalized','position',[.005 0.995],'verticalalignment','top','horizontalalignment','left','color','r') ;
        axis tight
        axis equal
        set(gca,'xtick',[],'ytick',[])
        set(gca,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
        set(gca,'ydir','reverse')
        box on
    for ii = 1:nFrames
        im.CData = IMG(:,:,:,ii) ;
        ttl.String = num2str(ii) ;
        drawnow
        %pause(0.01)
    end


    
%% PROCESS THE SEED

% Get the Seed
    Seed = hd.Seeds(seedNumber) ;

% Get Infos
    Elems = Seed.Triangles ;
    nNodesByElements = size(Elems,2) ;
    Nodes = Seed.Points ;
    nElems = size(Elems,1) ;
    nNodes = size(Nodes,1) ;

% PROCESSING ELEMENT GEOMETRY
    Xt = Nodes(:,1) ; Yt = Nodes(:,2) ; % Separate coordinates for indexing
    Xc = mean(Xt(Elems),2) ; Yc = mean(Yt(Elems),2) ; % Barycentric coordinates
    Areas = zeros(nElems,1) ; % Element areas
    gradT = zeros(nElems,2*nNodes,2,2) ; % Gradient matrix (element-wise)
    tri2nod = zeros(nNodes,nElems) ; % triangles to nodes transform
    wtbr = waitbar(0,'Processing elements') ;
    ticTime = tic ;
    for tri = 1:nElems
        % Points associated to the element
            Xtri = Nodes(Elems(tri,:),1) ;
            Ytri = Nodes(Elems(tri,:),2) ;
        % Sort nodes trigonometric orientation
            pol = angle(Xtri-Xc(tri)) + 1i*(Ytri-Yc(tri)) ;
            [~,ind] = sort(pol,'ascend') ;
            Elems(tri,:) = Elems(tri,ind) ;
            Xtri = Xtri(ind) ; Ytri = Ytri(ind) ;
        % Connectivity
            tri2nod(Elems(tri,:),tri) = 1 ;
        % Shape functions
            % Moments
                Areas(tri) = 1/2*det([ones(3,1) Xtri Ytri]) ;
                aa = cross(Xtri,Ytri) ;
                bb = circshift(Ytri,2) - circshift(Ytri,1) ;
                cc = circshift(Xtri,1) - circshift(Xtri,2) ;
        % Deformation gradient
            % Element-wise
                gradT(tri,Elems(tri,:),1,1) = 1/2/Areas(tri)*bb ; %dUx_dx
                gradT(tri,Elems(tri,:),1,2) = 1/2/Areas(tri)*cc ; %dUx_dy
                gradT(tri,nNodes+Elems(tri,:),2,1) = 1/2/Areas(tri)*bb ; %dUy_dx
                gradT(tri,nNodes+Elems(tri,:),2,2) = 1/2/Areas(tri)*cc ; %dUy_dy
        % waitbar
            if toc(ticTime)>0.02
                wtbr = waitbar(tri/nElems,wtbr) ;
                ticTime = tic ;
            end
    end
    wtbr = waitbar(0,wtbr,'Sparsing...') ;
    % TRANSFORM SOME OBJECTS TO SPARSE MATRICES
        gradT = {sparse(gradT(:,:,1,1)) sparse(gradT(:,:,1,2)) ; sparse(gradT(:,:,2,1)) sparse(gradT(:,:,2,2))} ;
        tri2nod = sparse(tri2nod) ;
        Valance = sum(tri2nod,2) ;
        invValance = sparse(1:nNodes,1:nNodes,1./Valance,nNodes,nNodes) ; % Inverse valance (for elemt averaging)

% BUILD SHAPE FUNCTIONS
    % For pixel-wise functions, sparse implementation
        totalArea = sum(abs(Areas)) ;
        nnzMap = round(totalArea)*10 ; % over-estimating the number of values
        iiiMAP = zeros(nnzMap,1) ; jjjMAP = zeros(nnzMap,1) ; vvvMAP = zeros(nnzMap,1) ; % Shape Functions
        iiiIN = zeros(nnzMap,1) ; jjjIN = zeros(nnzMap,1) ; vvvIN = zeros(nnzMap,1) ; % Inside elements
        lastIndiceMAP = 0 ;
        lastIndiceIN = 0 ;
    wtbr = waitbar(0,wtbr,'Building Shape Functions...') ;
    ticTime = tic ;
    for tri = 1:nElems
        % Points associated to the element
            Xtri = Nodes(Elems(tri,:),1) ;
            Ytri = Nodes(Elems(tri,:),2) ;
        % Pixels inside a rectangle around each element
            bbox = [floor(min(Xtri)) floor(min(Ytri)) ; ceil(max(Xtri)) ceil(max(Ytri))] ;
            IIe = II(bbox(1,2):bbox(2,2),bbox(1,1):bbox(2,1)) ;
            JJe = JJ(bbox(1,2):bbox(2,2),bbox(1,1):bbox(2,1)) ;
        % Pixels in the element
            IN = inpolygon(JJe,IIe,Xtri,Ytri) ;
            INDe = sub2ind([nI nJ],IIe(IN),JJe(IN)) ;
            nInd = length(INDe) ;
            iiiIN(lastIndiceIN+(1:nInd)) = INDe(:) ; % Pixel indices
            jjjIN(lastIndiceIN+(1:nInd)) = tri ; % Element
            vvvIN(lastIndiceIN+(1:nInd)) = 1 ; % Element
            lastIndiceIN = lastIndiceIN + nInd ;
        % Shape functions
            % Moments
                Areas(tri) = 1/2*det([ones(3,1) Xtri Ytri]) ;
                aa = cross(Xtri,Ytri) ;
                bb = circshift(Ytri,2) - circshift(Ytri,1) ;
                cc = circshift(Xtri,1) - circshift(Xtri,2) ;
            % Compute
                for pt = 1:3
                    iiiMAP(lastIndiceMAP+(1:nInd)) = INDe(:) ; % Pixel indices
                    jjjMAP(lastIndiceMAP+(1:nInd)) = Elems(tri,pt) ; % Node
                    vvvMAP(lastIndiceMAP+(1:nInd)) = (aa(pt) + bb(pt)*JJe(IN) + cc(pt)*IIe(IN))/2/Areas(tri) ; % Shape function value
                    lastIndiceMAP = lastIndiceMAP + nInd ;
                end
        % waitbar
            if toc(ticTime)>0.02
                wtbr = waitbar(tri/nElems,wtbr) ;
                ticTime = tic ;
            end
    end
    delete(wtbr)
    % INSIDE each Elements
        indValid = vvvIN~=0 ;
        INSIDE = sparse(iiiIN(indValid),jjjIN(indValid),vvvIN(indValid),nI*nJ,nElems) ;
        nInside = sum(INSIDE,1).' ;
        DOMAIN = logical(sum(INSIDE,2)) ;
    % SHAPE FUNCTIONS
        % Keep only the positive values and inferior to 1 
            indValid = vvvMAP~=0 ;
            %indValid = vvv>0 & vvv<=1 ;
        % Build the mapping matrix
            MAPPING = sparse(iiiMAP(indValid),jjjMAP(indValid),vvvMAP(indValid),nI*nJ,nNodes) ; % Shape functions (sparse representation)
    
    %DOMAIN = reshape(any(INSIDE,2),[nI nJ]) ;
    
%% COMPUTE CONSTANT OBJECTS

% Derivation Kernels
    N = 3 ;
    xx = (-N:N)' ;
    % Various kernel choices
        kern = [0 1 0]' ; dkern = [1 0 -1]'/2 ;
        %kern = cos(pi/2*xx/N).^2 ; dkern = -pi/N*sin(pi/2*xx/N).*cos(pi/2*xx/N) ;
        %sig = (N+1)/log(5*N) ; kern = exp(-xx.^2/sig^2) ; dkern =  -2*xx/sig^2.*exp(-xx.^2/sig^2) ;
    NORM = sum(sum(kern*kern')) ;
    Func = @(img)conv2(double(img),kern*kern','same')/NORM/imgClassRange ;
    dFunc_dx = @(img)conv2(double(img),kern*dkern','same')/NORM/imgClassRange ;
    dFunc_dy = @(img)conv2(double(img),dkern*kern','same')/NORM/imgClassRange ;

% Reference Image Processing
    % Convolution
        F = Func(img0) ;
    % Gradient of Reference image
        dF_dx = dFunc_dx(img0) ;
        dF_dy = dFunc_dy(img0) ;
        
% Projection of the Gradient in the nodal basis
    dF_da = [dF_dx(:).*MAPPING , dF_dy(:).*MAPPING] ;
    
% Hessian
    Hess = dF_da'*dF_da ;
    
    
%% STRAIN CRITERION

constraint = 'full' ; % 'full' or 'normal'
    
% Gradient matrix, [dUx_dx,dUy_dx,dUx_dy,dUy_dy]
    G = sparse(cat(1,gradT{1,1},gradT{2,1},gradT{1,2},gradT{2,2})) ;

% Linearized Green-Lagrane Epsilon = BB*a, epsilon 
    B = sparse(cat(1,gradT{1,1},gradT{2,2},0.5*(gradT{1,2}+gradT{2,1}))) ;
    BB = B'*B ; % To minimize the "strain energy"
    
% Minimisation of the strain GRADIENT in the dir. of the edge normal
    % Build edges associated to elements
        elmtEDGES = [Elems(:,1:2);Elems(:,2:3);Elems(:,[3,1])] ;
        elmtEDGES = sort(elmtEDGES,2) ;
        elmtEDGES = elmtEDGES(:,1)+1i*elmtEDGES(:,2) ; % Complex notation is more convenient
        elmtEDGES = reshape(elmtEDGES,[nElems 3]) ;
    % Edges of the mesh
        EDGES = unique(elmtEDGES) ;
        nEdges = length(EDGES) ;
    % Edge Connectivity
        edg2nod = sparse(nNodes,nEdges) ; % nodes linked to each edge
        tri2edg = sparse(nEdges,nElems) ; % elements linked to each edge
    % Edge constraint
        switch constraint
            case 'normal'
                E = sparse(2*nEdges,3*nElems) ; % projection on the normal
            case 'full'
                E = sparse(4*nEdges,3*nElems) ; % full gradient
        end
    % Process
        wtbr = waitbar(0,'Edge constraint...') ;
        for edg = 1:nEdges
            % Retrieve
                thisEdge = [real(EDGES(edg)) imag(EDGES(edg))] ;
                edg2nod(thisEdge,edg) = 1 ;
            % Find the attached triangle(s)
                [elmts,~] = find(abs(elmtEDGES-EDGES(edg))<eps) ;
                tri2edg(edg,elmts) = 1 ;
                if length(elmts)==1 ; continue ; end % the edge is naked
            % Edge normal
                tang = diff(Nodes(thisEdge,:),1,1) ; tang = tang.'/norm(tang) ;
                normal = [tang(1) tang(2) ; -tang(2) tang(1)]\[0;1] ;
                dist = diff([Xc(elmts(:))';Yc(elmts(:))'],1,2) ;
                dVect = [-1 1]/norm(dist) ;
            % Constrain
                switch constraint
                    case 'normal' % Projection of the STRAINS on the normal
                        % dExx.nx + dExy.ny = 0
                            E(edg,elmts) = normal(1)*dVect ;
                            E(edg,elmts+2*nElems) = normal(2)*dVect ;
                        % dExy.nx + dEyy.ny = 0
                            E(edg+nEdges,elmts+2*nElems) = normal(1)*dVect ;
                            E(edg+nEdges,elmts+1*nElems) = normal(2)*dVect ;
                    case 'full' % variation of the GRADIENT between elements
                        for comp = 1:4
                            E(edg+(comp-1)*nEdges,elmts+(comp-1)*nElems) = dVect ;
                        end
                end
            % Waitbar
                wtbr = waitbar(edg/nEdges,wtbr) ;
        end
        delete(wtbr) ;
    % Sparsify...
        E = sparse(E) ;
        edg2nod = sparse(edg2nod) ;
        tri2edg = sparse(tri2edg) ;
        
end % END OF INITIALIZATION
    
%% PERFORM DIC !

% PARAMETERS
    % Watch CPU 
        codeProfile = false ;
    % Plotting
        plotEachIteration = true ; % Plot at every iteration (without necessary pausing)
        pauseAtPlot = true ; % Pause at each iteration for debugging
        plotRate = Inf ; % Plot Refresh Frequency
    % Image Warping
        imWarpInterpOrder = 'cubic' ;
    % Image difference criterion
        diffCriterion = ... 'Diff' ... Simple difference
                        ... 'ZM_Diff' ... Zero-mean difference
                         'ZM_N_Diff' ... Normalized Zero-mean difference
                         ;
    % Geometry validation criteria
        cullOutOfFrame = true ; % Cull out of frame points
        localWEIGHT = INSIDE ; % MAPPING ; % For local averaging and difference image moments computations
        minCorrCoeff = .0 ; % Below this, elements are culled
        alwaysCheckCorrCoeff = false ; % Check corr. coeffs at each iteration or only after convergence ?
    % Regularization
        beta = 1*1e3 ; % Strain gradient penalisation coefficient
    % Convergence Criteria
        maxIt = 30 ; % Maximum number of Newton-Raphson iterations
        minNorm = 1e-4 ; % Maximum displacement of a node

% INIT FIGURE
    figure(figGlobalDIC) ;
    figGlobalDIC = clf(figGlobalDIC,'reset') ;
    figGlobalDIC.Tag = 'globalDICfigure' ;
        ax = [] ;
        ax(1) = mysubplot((nI<nJ)+1,(nI>=nJ)+1,1) ;
            im = imagesc(1:nJ,1:nI,Func(img0)) ;
            mesh = trisurf(Elems,Nodes(:,1),Nodes(:,2),Nodes(:,1)*0,'facecolor','none','edgecolor','r','linewidth',0.5,'edgealpha',0.5,'facealpha',0.5) ;
            markers = plot(NaN,NaN,'.b','markersize',15) ; % Deugging...
            colormap(ax(1),jet)
            set(ax(1),'Clipping','off') ;
        ax(2) = mysubplot((nI<nJ)+1,(nI>=nJ)+1,2) ;
            imRes = imagesc(1:nJ,1:nI,Func(img0)) ; 
            colormap(ax(2),gray)
            ttl = title('','interpreter','none','units','normalized','color',[1 0 0]*1.0,'position',[0 1],'horizontalalignment','left','verticalalignment','top') ;
        axis(ax,'tight')
        axis(ax,'equal')
        set(ax,'xtick',[],'ytick',[])
        set(ax,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
        set(ax,'ydir','reverse')
    stopBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','STOP'...
                        ,'units','normalized'...
                        ,'position',[0.01 0.01 .08 .05]...
                        ,'callback',@(src,evt)disp('Stop!')) ;
    if pauseAtPlot
        nextBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','NEXT'...
                        ,'units','normalized'...
                        ,'position',[0.1 0.01 .08 .05]...
                        ,'callback',@(src,evt)disp('Next!')) ;
        continueBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','CONTINUE'...
                        ,'units','normalized'...
                        ,'position',[0.19 0.01 .08 .05]...
                        ,'callback',@(src,evt)disp('Continue!')) ;
    end
    
% Other figure if needed to debug
    if 0
        if isempty(findobj(0,'tag','figDebug'))
            figDebug = figure('tag','figDebug') ;
        end
        figDebug = figure(findobj(0,'tag','figDebug')) ;
    end
        

% INITIALIZE
    % Nodes position
        Xn = ones([nNodes,2,nFrames])*NaN ;
        Xn(:,:,avgFrames) = repmat(Nodes,[1 1 length(avgFrames)]) ;
    % Displacements
        % Of Nodes
            Un = ones([nNodes,2,nFrames])*NaN ;
            Un(:,:,avgFrames) = 0 ;
        % Of Pixels
            Up = zeros([nI nJ 2]) ;
    % Reference Image
        img1 = F ;
        % Image Vector
            img1v = img1(:) ;
        % Moments
            % Integration weights
                sumWEIGHT = sum(localWEIGHT,1).' ;
            % Mean over elements
                meanImg1 = (localWEIGHT'*img1v)./sumWEIGHT(:) ;
            % Zero-local-mean on pixels
                img1m = img1v-localWEIGHT*meanImg1(:) ;
            % Norm over element
                normImg1 = sqrt(localWEIGHT'*(img1m(:).^2)) ;
            % Zero-local-mean-normalized images
                img1mz = img1m(:)./(localWEIGHT*normImg1) ;
    % Mask
        VALID = true(nNodes,1) ;
        validElems = true(nElems,1) ;
        validEdges = true(nEdges,1) ;
        nakedEdges = sum(tri2edg(:,validElems),2)<2  ;
    
% RUN !
if codeProfile ; profile on ; end
for ii = dicFrames
    % Import and display image
        img2 = Func(IMG(:,:,:,ii)) ;
        im.CData = repmat((img2-min(img2(:)))/range(img2(:)),[1 1 3]) ;
    % Init
        it = 0 ;
        outFlag = false ;
        Un(:,:,ii) = Un(:,:,ii-dicDir) ;
    % Add the previous "speed" as convergence help
        if abs(refFrame-ii)>=2
            Un(:,:,ii) = Un(:,:,ii) + (Un(:,:,ii-dicDir)-Un(:,:,ii-2*dicDir)) * 1.0 ;
        end
    % Newton-Raphson
        RMSE_0 = Inf ;
        lastPlotTime = tic ;
        while ~outFlag && ~stopBtn.Value
            % IMAGE WARPING
                % Compute the displacement at each pixel
                    Up = reshape(MAPPING(:,VALID)*Un(VALID,:,ii),[nI nJ 2]) ;
                % Warp the image (add 1i when the pixel is outside the frame)
                    img2w = interp2(JJ,II,img2,JJ+Up(:,:,1),II+Up(:,:,2),imWarpInterpOrder,1i) ;
            % VALID GEOMETRY
                deadPixels = imag(img2w(:))~=0 ;
                % Elements
                    validElems = validElems & sum(tri2nod(VALID,:),1)'==3 ; % triangles with still their three nodes valid
                    validElems = validElems & sum(INSIDE(~deadPixels,:),1)'~=0 ; % triangles with still some pixels in the image
                % Edges
                    validEdges = validEdges & sum(edg2nod(VALID,:),1)'==2 ; % edges with still their two ends points valid
                    nakedEdges = sum(tri2edg(:,validElems),2)<2 ; % edges with less than two elements valid
                % Nodes
                    VALID = VALID & sum(tri2nod(:,validElems),2)>0 ; % nodes must be linked to at least one valid element
                    VALID = VALID & sum(edg2nod(:,validEdges),2)>0 ; % nodes must be linked to at least one valid edge
                    nVALID = sum(VALID) ;
                    if nVALID==0 ; break ; end
            % VALID DIC Domain
                dicDomain = full(logical(sum(INSIDE(:,validElems),2))) ;
                dicDomain = dicDomain & ~deadPixels ;
            % CONVERT IMAGES TO VECTORS
                img2v = real(img2w(:).*dicDomain(:)) ;
            % IMAGE MOMENTS
                % Integration weights
                    switch size(localWEIGHT,2)
                        case nNodes
                            WEIGHT = localWEIGHT(:,VALID) ;
                            ww = 1./normImg1(VALID) ;
                        case nElems
                            WEIGHT = localWEIGHT(:,validElems) ;
                            ww = 1./bsxfun(@(x,y)x./y,tri2nod(VALID,validElems)*normImg1(validElems),sum(tri2nod(VALID,validElems),2)) ;
                    end
                    sumWEIGHT = sum(WEIGHT(dicDomain,:),1).' ;
                % Mean over elements
                    meanImg2 = (WEIGHT'*img2v)./sumWEIGHT(:) ;
                % Zero-local-mean on pixels
                    img2m = img2v-WEIGHT*meanImg2(:) ;
                % Norm over element
                    normImg2 = sqrt(WEIGHT'*(img2m(:).^2)) ;
                % Zero-local-mean-normalized images
                    img2mz = img2m(:)./(WEIGHT*normImg2) ;
            % IMAGE FUNCTIONAL
                switch diffCriterion
                    case 'Diff' % Simple difference
                        diffImg = img1v-img2v ;
                        weight = speye(2*nVALID) ;
                    case 'ZM_Diff' % Zero-mean difference
                        diffImg = img1m(:)-img2m(:) ;
                        weight = speye(2*nVALID) ;
                    case 'ZM_N_Diff' % Normalized Zero-mean difference
                        diffImg = img1mz(:)-img2mz(:) ;
                        weight = sparse(1:2*nVALID,1:2*nVALID,[ww;ww]) ;
                end
            % NEWTON-RAPHSON PROCEDURE
                % Compute the first RMSE derivative
                    dr_da = (dF_da(dicDomain,[VALID;VALID])'*diffImg(dicDomain)) ;
                % Contraint on the SECOND displacement gradient
                    switch constraint
                        case 'normal'
                            vEdg = repmat(~nakedEdges,[2 1]) ;
                            vEle = [validElems;validElems;validElems] ;
                            CONS = B(vEle,:)'*(E(vEdg,vEle)'*E(vEdg,vEle))*B(vEle,:) ;
                        case 'full'
                            vEdg = repmat(~nakedEdges,[4 1]) ;
                            vEle = [validElems;validElems;validElems;validElems] ;
                            CONS = G(vEle,:)'*(E(vEdg,vEle)'*E(vEdg,vEle))*G(vEle,:) ;
                    end
                    % Debugging...
                        %markers.XData = Xn(nodesOnNaked,1,ii); markers.YData = Xn(nodesOnNaked,2,ii);
                        %if any(~VALID)  disp('STOP!'); pause ; end
                % Updating DOFs, X*a=b
                    validDOF = [VALID;VALID] ;
                    X = ( ...
                            weight*Hess(validDOF,validDOF)*weight...
                            + beta*CONS(validDOF,validDOF)...
                        ) ; 
                    b = (...
                            weight*dr_da...
                            - beta*CONS(validDOF,validDOF)*[Un(VALID,1,ii);Un(VALID,2,ii)]...
                        ) ;
                    a = X\b ;
                % Residues
                    %residueImg = norm(Hess(validDOF,validDOF)*a-dr_da) ;
                    %residueCONS = norm(CONS(validDOF,validDOF)*a+CONS(validDOF,validDOF)*[Un(VALID,1,ii);Un(VALID,2,ii)]) ;
                % Displacement
                    Un(VALID,1,ii) = Un(VALID,1,ii) + a(1:nVALID) ;
                    Un(VALID,2,ii) = Un(VALID,2,ii) + a(nVALID+(1:sum(VALID))) ;
                % Positions
                    Xn(:,:,ii) = Nodes + Un(:,:,ii) ;
            % CONVERGENCE CITERIONS
                % Criterions
                    it = it+1 ;
                    %RMSE = norm(diffImg(DOMAIN))/norm(img1(DOMAIN)) ; 
                    normA = norm(a)/nNodes ; max(abs(a)) ;
                % Convergence criterion
                    if it>maxIt ; outFlag = true ; end
                    if normA<minNorm ; outFlag = true ; end
                    %if RMSE<1e-6 || abs((RMSE-RMSE_0)/RMSE) < 1e-4 ; outFlag = true ; end
                % Keep the error
                    %RMSE_0 = RMSE ;
            % POINTS/ELEMENTS VALIDATION/DELETION
                % OUT-OF-FRAME POINTS
                    if cullOutOfFrame
                        VALID = VALID & Xn(:,1,ii)<nJ+1 & Xn(:,1,ii)>0 & Xn(:,2,ii)<nI+1 & Xn(:,2,ii)>0 ;
                    end
                % DECORRELATED ELEMENTS
                    if minCorrCoeff>0 && (outFlag || alwaysCheckCorrCoeff)
                        corrCoeff = abs(WEIGHT'*(img1mz(:).*img2mz(:))) ;
                        %mesh.FaceVertexCData = zeros(nElems,1) ;
                        %mesh.FaceVertexCData(validElems) = corrCoeff(:) ;
                        %mesh.FaceColor = 'flat' ; caxis(ax(1),[0 1])
                        %;mesh.FaceAlpha = 1 ; colorbar(ax(1)) ; drawnow ;
                        switch size(localWEIGHT,2) 
                            case nElems % Correlation at the element level
                                %VALID = VALID & (tri2nod*corrCoeff(:))./sum(tri2nod(:,validElems),2) ;
                                validElems(validElems) = validElems(validElems) & corrCoeff(:)>minCorrCoeff ; 
                            case nNodes % Correlation at the node level
                                VALID(VALID) = VALID(VALID) & corrCoeff(:)>minCorrCoeff ;
                        end
                    end
                % SET THE NON-VALID NODES TO NAN
                    Xn(~VALID,:,ii) = NaN ;
                    Un(~VALID,:,ii) = NaN ;
            % DISPLAY
                if plotEachIteration || (outFlag && toc(lastPlotTime)>1/plotRate)
                    ttl.String = [num2str(ii),'(',num2str(it),')'] ;
                    imRes.CData = reshape(diffImg,[nI nJ]) ;
                    mesh.Vertices = Xn(:,:,ii) ;
                    %figure(figDebug) ; clf ; ind = (0:nJ-1)*nI+ceil(nI/2) ; plot(img1v(ind)) ; plot(img2v(ind)) ; 
                    drawnow ;
                    toc(lastPlotTime)
                    lastPlotTime = tic ;
                end
            % Pause execution ?
                if pauseAtPlot && ~continueBtn.Value
                    while ~stopBtn.Value && ~nextBtn.Value && ~continueBtn.Value
                        drawnow ;
                    end
                    nextBtn.Value = 0 ;
                    if continueBtn.Value
                        nextBtn.Visible = 'off' ;
                        continueBtn.Visible = 'off' ;
                    end
                end
        end
        % Out Criterions
            if stopBtn.Value; break ; end
            if nVALID==0; break ; end
end


% REVERSE THE DISPLACEMENTS IF NEEDED (backward mode)
    if dicDir<0
        firstValidFrame = sum(isnan(Xn(:,1,:)),3)+1 ;
        firstValidPosition = [...
                reshape(Xn(sub2ind(size(Xn),1:nNodes,ones(1,nNodes),firstValidFrame(:)')),[nNodes 1]) ...
                reshape(Xn(sub2ind(size(Xn),1:nNodes,2*ones(1,nNodes),firstValidFrame(:)')),[nNodes 1]) ...
                            ] ;
        Un = Xn-repmat(firstValidPosition,[1 1 nFrames]) ;
    end

% COMPUTE STRAINS
    % Valid elements as function of time
        VALID_t = reshape(~isnan(Un(:,1,:)),[nNodes nFrames]) ;
        validElems_t = false(nElems,nFrames) ;
        validElems_t(:,1:refFrame) = true ;
        for ii = 1:nFrames % no need to go to nImages if it has been stopped before
            validElems_t(:,ii) = sum(tri2nod(VALID_t(:,ii),:),1)'==3 ;
        end
    % Strains
        % DOFs as a function of time
            An = reshape(Un,[2*nNodes nFrames]) ; 
            An(isnan(An)) = 1i ; % Propagates NaNs with complex numbers
        % Gradient (at the element level)
            dUx_dx = gradT{1,1}*An + 1i*~validElems_t ;
            dUx_dy = gradT{1,2}*An + 1i*~validElems_t ;
            dUy_dx = gradT{2,1}*An + 1i*~validElems_t ;
            dUy_dy = gradT{2,2}*An + 1i*~validElems_t ;
        % Strains with NL terms
            Strains = [] ;
            Strains(:,1,:) = dUx_dx + 0.5*(dUx_dx.^2 + dUy_dx.^2) ;
            Strains(:,2,:) = dUy_dy + 0.5*(dUx_dy.^2 + dUy_dy.^2) ;
            Strains(:,3,:) = 0.5*(dUx_dy + dUy_dx + dUx_dx.*dUx_dy + dUy_dx.*dUy_dy) ;
            Strains(imag(Strains)~=0) = NaN ; % Re-set imaginary results to NaN
        % To the nodes level
            %meanOnElems = diag(1./sum(tri2nod(:,validElems),2))*tri2nod(:,validElems) ;
% SEND THE RESULT TO navDIC
        hd.Seeds(seedNumber).Strains = interpn(...
                                                repmat((1:size(Strains,1))',[1 3 nFrames]),...
                                                repmat(1:3,[size(Strains,1) 1 nFrames]),...
                                                repmat(reshape(frames,[1 1 nFrames]),[size(Strains,1) 3 1]),...
                                                Strains,...
                                                repmat((1:size(Strains,1))',[1 3 hd.nFrames]),...
                                                repmat(1:3,[size(Strains,1) 1 hd.nFrames]),...
                                                repmat(reshape(navDICFrames,[1 1 hd.nFrames]),[size(Strains,1) 3 1]),...
                                            'linear',NaN) ;
        hd.Seeds(seedNumber).MovingPoints = interpn(...
                                                repmat((1:nNodes)',[1 2 nFrames]),...
                                                repmat(1:2,[nNodes 1 nFrames]),...
                                                repmat(reshape(frames,[1 1 nFrames]),[nNodes 2 1]),...
                                                Xn,...
                                                repmat((1:nNodes)',[1 2 hd.nFrames]),...
                                                repmat(1:2,[nNodes 1 hd.nFrames]),...
                                                repmat(reshape(navDICFrames,[1 1 hd.nFrames]),[nNodes 2 1]),...
                                            'linear',NaN) ;
        hd.Seeds(seedNumber).Displacements = interpn(...
                                                repmat((1:nNodes)',[1 2 nFrames]),...
                                                repmat(1:2,[nNodes 1 nFrames]),...
                                                repmat(reshape(frames,[1 1 nFrames]),[nNodes 2 1]),...
                                                Un,...
                                                repmat((1:nNodes)',[1 2 hd.nFrames]),...
                                                repmat(1:2,[nNodes 1 hd.nFrames]),...
                                                repmat(reshape(navDICFrames,[1 1 hd.nFrames]),[nNodes 2 1]),...
                                            'linear',NaN) ;
    
    
if codeProfile; profile viewer ; profile off ; end

return ;









%% REVERSE THE DISPLACEMENTS (3D PRINTING IMAGES)

% New reference state: just before the 'disparition' of the point
lastValidFrame = sum(~isnan(Xn(:,1,:)),3) ;
lastValidPosition = [...
        reshape(Xn(sub2ind(size(Xn),1:nNodes,ones(1,nNodes),lastValidFrame(:)')),[nNodes 1]) ...
        reshape(Xn(sub2ind(size(Xn),1:nNodes,2*ones(1,nNodes),lastValidFrame(:)')),[nNodes 1]) ...
                    ] ;
revUn = Xn-repmat(lastValidPosition,[1 1 nImages]) ;

% Cull points with a maximum disp. Criterion
    maxDispCrit = 20 ; Inf ;
    maxDisp = sqrt(max(sum(revUn.^2,2),[],3)) ;
    revUn(maxDisp>maxDispCrit,:,:) = NaN ;

% SEND THE RESULT TO navDIC

% Seed number
    seedNumber = 1 ;

% Displacements/Positions
    hd.Seeds(seedNumber).MovingPoints = Xn ;
    hd.Seeds(seedNumber).Displacements = revUn ;
    
% Strains
    An = reshape(revUn,[2*nNodes nImages]) ; 
    An(isnan(An)) = 1i ; % Propagates false positives
    dUx_dx = invValance*tri2nod*gradT{1,1}*An ;
    dUx_dy = invValance*tri2nod*gradT{1,2}*An ;
    dUy_dx = invValance*tri2nod*gradT{2,1}*An ;
    dUy_dy = invValance*tri2nod*gradT{2,2}*An ;
    Strains = [] ;
    Strains(:,1,:) = dUx_dx + 0.5*(dUx_dx.^2 + dUy_dx.^2) ;
    Strains(:,2,:) = dUy_dy + 0.5*(dUx_dy.^2 + dUy_dy.^2) ;
    Strains(:,3,:) = 0.5*(dUx_dy + dUy_dx + dUx_dx.*dUx_dy + dUy_dx.*dUy_dy) ;
    Strains(imag(Strains)~=0) = NaN ; % Re-set to NaN
    hd.Seeds(seedNumber).Strains = Strains ;
    
    
if codeProfile; profile viewer ; profile off ; end

    
%% PREVIOUS EXPORT FUNCTIONS
    
% Strains in each triangle
    Et = zeros(nElems,3,nImages) ; % Triangle level
    Btot = zeros(nElems,2*nNodes,3) ;
    Bp = zeros(nNodes,2*nNodes,3) ;
    Ep = zeros(nNodes,3,nImages) ; % Point level
        for tri = 1:nElems
            pts = Elems(tri,:) ;
            Xt = Nodes(pts,:) ;
            A = 1/2 * det([ones(3,1) Xt]) ;
            B = 1/2/A*[Xt(2,2)-Xt(3,2) 0 Xt(3,2)-Xt(1,2) 0 Xt(1,2)-Xt(2,2) 0;...
                       0 Xt(3,1)-Xt(2,1) 0 Xt(1,1)-Xt(3,1) 0 Xt(2,1)-Xt(1,1);...
                       Xt(3,1)-Xt(2,1) 0 Xt(1,1)-Xt(3,1) 0 Xt(2,1)-Xt(1,1) 0;...
                       0 Xt(2,2)-Xt(3,2) 0 Xt(3,2)-Xt(1,2) 0 Xt(1,2)-Xt(2,2)] ;
            for ii = 1:nImages
                Ut = Un(pts,:,ii).' ;
                epsilon = B*Ut(:) ;
                epsilon = [epsilon(1) epsilon(3) ; epsilon(4) epsilon(2)] ;
                %epsilon = 0.5*(epsilon + epsilon' + epsilon' * epsilon) ;
                epsilon = [epsilon(1,1) epsilon(2,2) epsilon(1,2)+epsilon(2,1)] ;
                Et(tri,:,ii) = epsilon(1:3) ;
            end
            dofs = (2*(pts(:)-1)+(1:2))' ;
            Btot(tri,dofs,:) = (B(1:3,:) + B([1 2 4],:))' ;
        end
    % Point level
        for pt = 1:nNodes
            [tris,~] = find(Elems==pt) ;
            Ep(pt,:,:) = mean(Et(tris,:,:),1) ;
            Bp(pt,:,:) = mean(Btot(tris,:,:),1) ;
        end
    % Save
        Seed.Strains = Ep ;
        hd.Seeds(seedNumber) = Seed ;



    
%% SVD ON THE MEASURED DISPLACEMENT

Usvd = reshape(Un(:,:,1:10),nNodes*2,[]) ;
nSamples = size(Usvd,2) ;

Css = Usvd*Usvd' ;
    
[u,s,v] = svd(Css) ;
Ueig = u(:,1:nSamples) ;
Veig = v(:,1:nSamples) ;
S = diag(s(1:nSamples,1:nSamples)) 

clf
%plot(1-cumsum(diag(s))/sum(diag(s))),set(gca,'yscale','log')
plot(S),set(gca,'yscale','log')
    
clf
colormap jet
nModes = 5 ;
amp = norm([nI nJ]) ;
for m = 1:nModes
    mysubplot(ceil(nModes/2),2,m) ;
    trisurf(Elems,Nodes(:,1),Nodes(:,2),0*Nodes(:,1),'facecolor','none','edgecolor','k')
    tr = trisurf(Elems,Nodes(:,1),Nodes(:,2),0*Nodes(:,1),'facecolor','interp','edgecolor','none','facealpha',.5) ;
    tr.Vertices = tr.Vertices + amp*[Ueig(1:nNodes,m) Ueig((1:nNodes)+nNodes,m) tr.Vertices(:,3)*0] ;
    %tr.CData = sqrt(Ueig(1:nNodes,m).^2 + Ueig((1:nNodes)+nNodes,m).^2) ;
    tr.CData = Bp(:,:,3)*reshape(reshape(Ueig(:,m),[],2)',[],1) ;
    axis off
    axis equal
    set(gca,'xtick',[],'ytick',[])
    axis tight
end


%% EFFECT OF THE AVERAGING OF THE FIRST IMAGES

clf reset ;
axes('position',[0 0 1 1])
    im = imagesc(1:nJ,1:nI,img0) ; colormap(gray)
    ttl = title('','interpreter','none') ;
    axis tight
    axis equal
    set(gca,'xtick',[],'ytick',[])
    set(gca,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
    box on
    %caxis([0 1])

img = img0 ;
for ii = 1:85%nImages
    img = IMG{ii};%(img*ii + IMG{ii})/(ii+1) ;
    im.CData = img ;
    ttl.String = num2str(ii) ;
    drawnow ;
end



    
%% TEST THE GRADIENT OPERATOR
    gradTh = [1 0 ; 0 0]*1e-1
    dx = Nodes(:,1)-mean(Nodes(:,1)) ; dy = Nodes(:,2)-mean(Nodes(:,2)) ;
    U1 = dx.^2/1000*gradTh(1,1) + dy*gradTh(1,2) ; U2 = dx*gradTh(2,1) + dy*gradTh(2,2) ;
    UU = An(:,12) ;[U1 ; U2] ;
    dUx_dx = gradT(:,:,1,1)*UU ;
    dUx_dy = gradT(:,:,1,2)*UU ;
    dUy_dx = gradT(:,:,2,1)*UU ;
    dUy_dy = gradT(:,:,2,2)*UU ;
    gradResTris = [norm(dUx_dx-gradTh(1,1)) norm(dUx_dy-gradTh(1,2)) ; norm(dUy_dx-gradTh(2,1)) norm(dUy_dy-gradTh(2,2))] %
    gradResNodes = [norm(diag(1./Valance)*tri2nod*dUx_dx-gradTh(1,1)) norm(diag(1./Valance)*tri2nod*dUx_dy-gradTh(1,2)) ; norm(diag(1./Valance)*tri2nod*dUy_dx-gradTh(2,1)) norm(diag(1./Valance)*tri2nod*dUy_dy-gradTh(2,2))] %
    
    
    clf reset
    %trimesh(Elems,Nodes(:,1),Nodes(:,2),0*U2,'edgecolor','k','facecolor','none') ;
    if 0 % Facet strain
        trisurf(Elems,Nodes(:,1)+U1,Nodes(:,2)+U2,0*U2,gradT(:,:,2,2)*UU,'facealpha',0.5,'facecolor','flat') ;
    else % Point strain
        trisurf(Elems,Nodes(:,1)+U1,Nodes(:,2)+U2,0*U2,diag(1./Valance)*tri2nod*gradT(:,:,2,2)*UU,'facealpha',0.5,'facecolor','interp') ;
    end
    axis equal
    axis tight
    colorbar
    
    
    

    
%% VERIFY MAPPING FUNCTIONS
    if 1
        clf reset ;
        axes('position',[0 0 1 1])
            im = imagesc(1:nJ,1:nI,img0) ; colormap(flip(gray))
            ttl = title('','interpreter','none') ;
            trisurf(Elems,Nodes(:,1),Nodes(:,2),0*Nodes(:,2),'facecolor','none','edgecolor','k')
            axis tight
            axis equal
            %set(gca,'xtick',[],'ytick',[])
            set(gca,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
            set(gca,'ydir','reverse')
            box on
            caxis([0 1])
        for nod = 1:nNodes
            im.CData = reshape(MAPPING(:,nod),[nI nJ]) ;
            ttl.String = num2str(nod) ;
            drawnow
            %pause(0.01)
        end
    end

    