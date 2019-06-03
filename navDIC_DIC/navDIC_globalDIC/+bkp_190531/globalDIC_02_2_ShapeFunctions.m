
    
% BUILD SHAPE FUNCTIONS
    % For pixel-wise functions, sparse implementation
        totalArea = sum(abs(Areas)) ;
        nnzMap = round(totalArea)*10 ; % over-estimating the number of values
        iiiMAP = zeros(nnzMap,1) ; jjjMAP = zeros(nnzMap,1) ; vvvMAP = zeros(nnzMap,1) ; % Shape Functions
        iiiIN = zeros(nnzMap,1) ; jjjIN = zeros(nnzMap,1) ; vvvIN = false(nnzMap,1) ; % Inside elements
        lastIndiceMAP = 0 ;
        lastIndiceIN = 0 ;
    wtbr = waitbar(0,'Building Shape Functions...') ;
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
            vvvIN(lastIndiceIN+(1:nInd)) = true ; % Element
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
    wtbr = waitbar(0,wtbr,'Sparsing...') ;
    % INSIDE each Elements
        INSIDE = sparse(iiiIN(vvvIN),jjjIN(vvvIN),true,nI*nJ,nElems) ;
        nInside = sum(INSIDE,1).' ;
        DOMAIN = logical(sum(INSIDE,2)) ;
    % SHAPE FUNCTIONS
        % Keep only the positive values and inferior to 1 
            indValid = vvvMAP~=0 ;
            %indValid = vvv>0 & vvv<=1 ;
        % Build the mapping matrix
            MAPPING = sparse(iiiMAP(indValid),jjjMAP(indValid),vvvMAP(indValid),nI*nJ,nNodes) ; % Shape functions (sparse representation)

    delete(wtbr)