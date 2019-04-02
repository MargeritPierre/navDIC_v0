clc
global hd
clearvars -except hd

% TIME THE CODE ?
codeProfile = false ;
profile off
if codeProfile; profile on ; end

% RETRIEVE THE DATA
    IMG = hd.Images{1}(:,:,:,1:1:end) ;
 
% PROCESS IMAGES /!\ STILL NO SUPPORT FOR MULTICOLOR IMAGES YET
    nImages = size(IMG,4) ;
    img0 = IMG(:,:,:,1) ;
    [nI,nJ] = size(img0) ;
    [JJ,II] = meshgrid(1:nJ,1:nI) ;
    %IMG = simple(IMG) ;
    for ii = 1:nImages
        %IMG(:,:,:,ii) = IMG(:,:,:,ii)/norm(IMG{ii}(:)) ;
        %IMG(:,:,:,ii) = IMG{ii}-mean(IMG{ii}(:)) ;
    end
    img0 = IMG(:,:,:,1) ;
    
% TIME MEANING
    timeMeanLength = 0 ;
    if timeMeanLength>0
        meanTime = 2*timeMeanLength+1 ;
        TimeKernel = ones(meanTime,1)/meanTime ;
        IMG0 = IMG ;
        IMG = IMG*0 ;
        wtbr = waitbar(0,'Time Averaging') ;
        for ii = 1:nImages
            for tt = -timeMeanLength:timeMeanLength
                ind = max(1,min(nImages,ii+tt)) ;
                IMG(:,:,:,ii) = IMG(:,:,:,ii) + IMG0(:,:,:,ind)*TimeKernel(tt+timeMeanLength+1) ;
            end
            wtbr = waitbar(ii/nImages,wtbr,['Time Averaging (',num2str(ii),'/',num2str(nImages),')']) ;
        end
        delete(wtbr) ;
        clear IMG0
    end

% REFERENCE IMAGE
    firstImage = 1 ;
    if 1 % AVERAGING THE FIRST IMAGES
        img0 = IMG(:,:,:,1)*0 ;
        for ii = 1:firstImage
            img0 = img0*((ii-1)/ii) + IMG(:,:,:,ii)/ii ;
        end
    end
    
%% PLAY VIDEO OF PROCESSED IMAGES

    clf reset ;
    axes('position',[0 0 1 1])
        im = imagesc(1:nJ,1:nI,img0) ; colormap(gray)
        ttl = title('','interpreter','none','units','normalized','position',[.005 0.995],'verticalalignment','top','horizontalalignment','left','color','r') ;
        axis tight
        axis equal
        set(gca,'xtick',[],'ytick',[])
        set(gca,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
        set(gca,'ydir','reverse')
        box on
    for ii = 1:nImages
        im.CData = IMG(:,:,:,ii) ;
        ttl.String = num2str(ii) ;
        drawnow
        %pause(0.01)
    end


    
%% PROCESS THE SEED

% Get the Seed
    Seed = hd.Seeds(1) ;

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
    Func = @(img)conv2(img,kern*kern','same')/NORM ;
    dFunc_dx = @(img)conv2(img,kern*dkern','same')/NORM ;
    dFunc_dy = @(img)conv2(img,dkern*kern','same')/NORM ;

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
            case "normal"
                E = sparse(2*nEdges,3*nElems) ; % projection on the normal
            case "full"
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
                dVect = [-1 1]*sum(normal.*dist) ;
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
    
%% PERFORM DIC !

% INIT FIGURE
    plotEachIteration = true ;
    fig = findobj(groot,'tag','globalDICfigure') ;
        if isempty(fig)
            fig = figure ;
        end
    fig = clf(fig,'reset') ;
    fig.Tag = 'globalDICfigure' ;
        ax = [] ;
        ax(1) = mysubplot((nI<nJ)+1,(nI>=nJ)+1,1) ;
            im = imagesc(1:nJ,1:nI,Func(img0)) ;
            mesh = trisurf(Elems,Nodes(:,1),Nodes(:,2),Nodes(:,1)*0,'facecolor','none','edgecolor','r','linewidth',0.5,'edgealpha',0.5,'facealpha',0.5) ;
            markers = plot(NaN,NaN,'.b','markersize',15) ; % Deugging...
            colormap(ax(1),jet)
        ax(2) = mysubplot((nI<nJ)+1,(nI>=nJ)+1,2) ;
            imRes = imagesc(1:nJ,1:nI,Func(img0)) ; 
            colormap(ax(2),gray)
            ttl = title('','interpreter','none','units','normalized','color',[1 0 0]*1.0,'position',[0 1],'horizontalalignment','left','verticalalignment','top') ;
        axis(ax,'tight')
        axis(ax,'equal')
        set(ax,'xtick',[],'ytick',[])
        set(ax,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
        set(ax,'ydir','reverse')
    stopBtn = uicontrol(fig,'style','togglebutton'...
                        ,'string','STOP'...
                        ,'units','normalized'...
                        ,'position',[0.01 0.01 .08 .05]...
                        ,'callback',@(src,evt)disp('click!')) ;


% INITIALIZE
    % Nodes position
        Xn = ones([nNodes,2,nImages])*NaN ;
        Xn(:,:,1:firstImage) = repmat(Nodes,[1 1 firstImage]) ;
    % Displacements
        % Of Nodes
            Un = ones([nNodes,2,nImages])*NaN ;
            Un(:,:,1:firstImage) = 0 ;
        % Of Pixels
            Up = zeros([nI nJ 2]) ;
    % Images
        img1 = F ;
    % Mask
        VALID = true(nNodes,1) ;
        validElems = true(nElems,1) ;
        validEdges = true(nEdges,1) ;
        nakedEdges = sum(tri2edg(:,validElems),2)<2  ;
    
% RUN !
for ii = firstImage+1:nImages
    % Import and display image
        img2 = Func(IMG(:,:,:,ii)) ;
        im.CData = repmat((img2-min(img2(:)))/range(img2(:)),[1 1 3]) ;
    % Init
        it = 0 ;
        outFlag = false ;
        Un(:,:,ii) = Un(:,:,ii-1) ;
    % Add the previous "speed" as convergence help
        if ii>2
            Un(:,:,ii) = Un(:,:,ii) + (Un(:,:,ii-1)-Un(:,:,ii-2)) * 1.0 ;
        end
    % Newton-Raphson
        RMSE_0 = Inf ;
        lastPlotTime = tic ;
        while ~outFlag && ~stopBtn.Value
            % VALID GEOMETRY
                % Elements
                    validElems = validElems & sum(tri2nod(VALID,:),1)'==3 ; % triangles with still their three nodes valid
                % Edges
                    validEdges = validEdges & sum(edg2nod(VALID,:),1)'==2 ; % edges with still their two ends points valid
                    nakedEdges = sum(tri2edg(:,validElems),2)<2 ; % edges with less than two elements valid
                % Nodes
                    VALID = VALID & sum(tri2nod(:,validElems),2)>0 ; % nodes must be linked to at least one valid element
                    VALID = VALID & sum(edg2nod(:,validEdges),2)>0 ; % nodes must be linked to at least one valid edge
                    nVALID = sum(VALID) ;
                    if nVALID==0 ; break ; end
                % DIC Domain
                    dicDomain = logical(sum(INSIDE(:,validElems),2)) ;
            % IMAGE WARPING
                % Compute the displacement at each pixel
                    Up = reshape(MAPPING(:,VALID)*Un(VALID,:,ii),[nI nJ 2]) ;
                % Warp the image
                    img2w = interp2(JJ,II,img2,JJ+Up(:,:,1),II+Up(:,:,2),'cubic',0) ;
            % IMAGE MOMENTS
                WEIGHT = INSIDE ; % MAPPING ; %
                sumWEIGHT = sum(WEIGHT,1).' ;
                % Mean over elements
                    meanImg1  = (WEIGHT'*img1(:))./sumWEIGHT(:) ;
                    meanImg2w  = (WEIGHT'*img2w(:))./sumWEIGHT(:) ;
                % Zero-local-mean on pixels
                    img1m = img1(:)-WEIGHT*meanImg1(:) ;
                    img2wm = img2w(:)-WEIGHT*meanImg2w(:) ;
                % Norm over element
                    normImg1 = sqrt(WEIGHT'*(img1m(:).^2)) ;
                    normImg2w = sqrt(WEIGHT'*(img2wm(:).^2)) ;
                % Zero-local-mean-normalized images
                    img1mz = img1m(:)./(WEIGHT*normImg1) ;
                    img2wmz = img2wm(:)./(WEIGHT*normImg2w) ;
            % IMAGE FUNCTIONAL
                diffImg = img1(:)-img2w(:) ; % Simple difference
                %diffImg = img1m(:)-img2wm(:) ; % Zero-mean difference
                %diffImg = img1mz(:)-img2wmz(:) ; % Normalized Zero-mean difference
            % NEWTON-RAPHSON PROCEDURE
                % Compute the first RMSE derivative
                    dr_da = (diffImg(dicDomain)'*dF_da(dicDomain,[VALID;VALID]))' ;
                % Contraint on the SECOND displacement gradient
                    beta = 0*1e-5 ; % penalisation coefficient
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
                % Updating DOFs
                    validDOF = [VALID;VALID] ;
                    a = ( ...
                            Hess(validDOF,validDOF)...
                            + beta*CONS(validDOF,validDOF)...
                        )\(...
                            dr_da...
                            - beta*CONS(validDOF,validDOF)*[Un(VALID,1,ii);Un(VALID,2,ii)]...
                        ) ;
                % Residues
                    %residueImg = norm(Hess(validDOF,validDOF)*a-dr_da) ;
                    %residueCONS = norm(CONS(validDOF,validDOF)*a+CONS(validDOF,validDOF)*[Un(VALID,1,ii);Un(VALID,2,ii)]) ;
                % Displacement
                    Un(VALID,1,ii) = Un(VALID,1,ii) + a(1:nVALID) ;
                    Un(VALID,2,ii) = Un(VALID,2,ii) + a(nVALID+(1:sum(VALID))) ;
                % Positions
                    Xn(:,:,ii) = Nodes + Un(:,:,ii) ;
            % CORRELATION COEFFICIENT
                % CULL OUT-OF-FRAME POINTS
                    VALID = VALID & Xn(:,1,ii)<nJ+1 & Xn(:,1,ii)>0 & Xn(:,2,ii)<nI+1 & Xn(:,2,ii)>0 ;
                % Decorrelated elements
                    minCorrCoeff = 0.1 ;
                    if minCorrCoeff>0
                        corrCoeff = abs(WEIGHT'*(img1mz(:).*img2wmz(:))) ;
                        switch size(WEIGHT,2) 
                            case nElems % Correlation at the element level
                                %VALID = VALID & (tri2nod*corrCoeff(:))./sum(tri2nod(:,validElems),2) ;
                                validElems = validElems & corrCoeff(:)>minCorrCoeff ; 
                            case nNodes % Correlation at the node level
                                VALID = VALID & corrCoeff(:)>minCorrCoeff ;
                        end
                    end
                % Set the non-valid values to NaN
                    Xn(~VALID,:,ii) = NaN ;
                    Un(~VALID,:,ii) = NaN ;
            % CONVERGENCE CITERIONS
                % Criterions
                    it = it+1 ;
                    %RMSE = norm(diffImg(DOMAIN))/norm(img1(DOMAIN)) ; 
                    normA = norm(a)/nNodes ; max(abs(a)) ;
                % Convergence criterion
                    if it>20 ; outFlag = true ; end
                    if normA<1e-4 ; outFlag = true ; end
                    %if RMSE<1e-6 || abs((RMSE-RMSE_0)/RMSE) < 1e-4 ; outFlag = true ; end
                % Keep the error
                    %RMSE_0 = RMSE ;
            % DISPLAY
                if plotEachIteration || (outFlag && toc(lastPlotTime)>0.05)
                    ttl.String = [num2str(ii),'(',num2str(it),')'] ;
                    imRes.CData = reshape(diffImg,[nI nJ]) ;
                    mesh.Vertices = Xn(:,:,ii) ;
                    if 0 % Show Strains
                        % Compute the strains
                            gradU = reshape(G(:,validDOF)*[Un(VALID,1,ii);Un(VALID,2,ii)],[nElems 2 2]) ;
                            epsU = 0.5*(gradU + permute(gradU,[1 3 2]) + permute(sum(permute(gradU,[1 3 2]).*permute(gradU,[1 4 2 3]),3),[1 2 4 3])) ;
                        mesh.FaceVertexCData = invValance*tri2nod*epsU(:,1,2) ; mesh.FaceColor = 'interp' ; %colorbar(ax(1))
                    end
                    %mesh.FaceVertexCData = corrCoeff(:) ; mesh.FaceColor = 'flat' ; caxis(ax(1),[0 1]) ;mesh.FaceAlpha = 1 ; colorbar(ax(1))
                    drawnow ;
                    lastPlotTime = tic ;
                end
            % Pause execution ?
                %pause
        end
        % Out Criterions
            if stopBtn.Value; break ; end
            if nVALID==0; break ; end
end


% SEND THE RESULT TO navDIC

% Seed number
    seedNumber = 1 ;

% Displacements/Positions
    hd.Seeds(seedNumber).MovingPoints = Xn ;
    hd.Seeds(seedNumber).Displacements = Un ;
    
% Strains
    An = reshape(Un,[2*nNodes nImages]) ; 
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
        hd.Seeds(1) = Seed ;



    
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

    