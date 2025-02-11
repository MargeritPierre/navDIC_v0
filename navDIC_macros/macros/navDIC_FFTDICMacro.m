classdef navDIC_FFTDICMacro < navDIC_DICMacro
%NAVDIC_FFTDICMacro Implements local DIC via FFT

properties % see fftdispmethod script for more details
    CorrSize = 0 ; % Rectangular correlation window (imaget size) see this.autoCorrSize
    Margin = 1/4 ; % Margin to truncate borders (relative to window size)
    MaxImagetShift = 1/2 ; % Maximum allowed imagette shift (close to image borders) (relative to window size)
    MaxDispPerIteration = 1/2 ; % Maximum allowed displacement per iteration (relative to window size)
    MaxDispAtConvergence = 1.1 ; % continue itrating if update is more than n pixels
    MaxIterations = 10 ; % maximum number of iterations
    Method = 'COR' ; % type of correlation
    Windowing = true ; % apply a blackman windowing
end

methods
    function this = navDIC_FFTDICMacro()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setupUI(this,hd)
    % Setup the macro (change parameters, etc)
        hd = setupUI@navDIC_DICMacro(this,hd) ; % Set parameters common to DIC macros
        if isempty(this.Seed) ; return ; end
        this.CorrSize = this.autoCorrSize() ;
    % Set globalDIC parameters
        defInputs = { ...
                         'Correlation window size [nX nY] 0:auto, -1:full' , mat2str(this.CorrSize) ...
                         ; 'FFT spectrum trunctation (relative)' , num2str(this.Margin) ...
                         ; 'Maximum imaget shift (relative)' , num2str(this.MaxImagetShift) ...
                         ; 'Maximum displacement per iteration (relative)' , num2str(this.MaxDispPerIteration) ...
                         ; 'Maximum displacement at convergence (pixels)' , num2str(this.MaxDispAtConvergence) ...
                         ; 'Maximum number of iterations' , num2str(this.MaxIterations) ...
                         ; 'Correlation method [COR/LS/SVD]' , this.Method ...
                         ; 'Apply windowing' , num2str(this.Windowing) ...
                    } ;
        out = inputdlg(defInputs(:,1),'FFT DIC parameters',1,defInputs(:,2)) ;
        if isempty(out) ; return ; end
        this.CorrSize = str2num(out{1}) ; 
        this.Margin = str2num(out{2}) ;
        this.MaxImagetShift = str2num(out{3}) ; 
        this.MaxDispPerIteration = str2num(out{4}) ; 
        this.MaxDispAtConvergence = str2num(out{5}) ; 
        this.MaxIterations = str2num(out{6}) ;
        this.Method = out{7} ; 
        this.Windowing = str2num(out{8}) ;
    % Setup the macro
        hd = this.setup(hd) ;
    end
    
    function hd = setup(this,hd)
    % Prepare the DIC data
        hd = setup@navDIC_DICMacro(this,hd) ; 
        this.setupDIC(hd) ;
    end
end


%% DIC PROCEDURES
methods
    function sz = autoCorrSize(this)
    % Automatic choice of correlation window sizes
        res = size(this.RefImgs{end},[2 1]) ;
        switch sum(abs(this.CorrSize))
            case 0 % automatic choice based on seed elements
                switch size(this.Seed.Elems,2)
                    case 1 % Points, take the frame size divided by he number of points
                        nPts = size(this.Seed.Points,1) ;
                        if nPts==1 
                            sz = res-2 ;
                        else
                            nPix = prod(res) ;
                            nPixByPt = nPix./nPts ;
                            sz = floor(sqrt(nPixByPt)).*[1 1] ;
                        end
                    case 2 % Edges, take the mean edge length
                        edgPts = reshape(this.Seed.Points(this.Seed.Elems,:),[],2,2) ;
                        Le = sqrt(sum(diff(edgPts,1,2).^2,3)) ;
                        sz = floor(median(Le)).*[1 1] ;
                    otherwise % tiangles or quads, use the mean element area
                        [~,A] = this.Seed.jacobian([.5 .5]) ;
                        sz = floor(sqrt(median(A))).*[1 1] ;
                end
            case -1 % Full image size
                sz = res - 2 ;
            otherwise % user-defined choice [nX nY]
                sz = round(this.CorrSize) ;
        end
    end
    
    function setupDIC(this,hd)
    % Setup all variables for the DIC...
        this.CorrSize = this.autoCorrSize() ;
    end
    
    function X = updateDIC(this,X,imgs,refImgs)
    % Update the configuration X using DIC performed on imgs
    % For now, only working on one gray-scale image
        G = refImgs{end}(:,:,1) ; % reference <MONOCHROME SINGLE CAMERA!>
        g = imgs{end}(:,:,1) ; % current <MONOCHROME SINGLE CAMERA!>
        X0 = this.Seed.Points ;
        params = struct(...
                         'dir' , this.DispComp ...
                       , 'CorrSize' , this.CorrSize ...
                       , 'm' , round(this.Margin.*this.CorrSize) ...
                       , 'maxImagetShift' , round(this.MaxImagetShift.*this.CorrSize) ...
                       , 'uMax' , round(this.MaxDispPerIteration.*this.CorrSize) ...
                       , 'iterateIfMaxDispHigherThan' , this.MaxDispAtConvergence ...
                       , 'maxIt' , this.MaxIterations ...
                       , 'FIT' , this.Method ...
                       , 'windowing' , this.Windowing ...
                        ) ;
        X = fftDispMethod(X,X0,g,G,params) ;
    end
    
    function I = transformImage(this,I,x,X)
    % Transform the image(s) I from configuration X to x
        x = reshape(this.N2V*x(:),[],2) ; 
        for ii = 1:numel(I)
            I{ii}(this.ROI) = this.bilinearInterp(I{ii},this.MAP*x) ;
        end
    end
end

%% MESH FUNCTIONS
methods
    function [Vertices,Tris,N2V] = simplexMesh(this,Nodes,Elems)
    % Build an augmented simplex mesh from an heterogeneous mesh
    % Elems can be edges, triangles, quads, or polygonal
        if nargin<2 ; Nodes = this.Seed.Points ; end % reference nodes
        if nargin<3 ; Elems = this.Seed.Elems ; end % seed elements
        
    % Infos
        nNodes = size(Nodes,1) ;
        
    % Original nodes
        Vertices = Nodes ;
        N2V = speye(2*nNodes) ;
        
    % Original triangle elements
        isTri = sum(Elems>0,2)==3 ; % not nan or with valid index (>0)
        if any(isTri) ; Tris = Elems(isTri,1:3) ;
        else ; Tris = zeros(0,3) ; end
        
    % Other elements
        Elems = Elems(~isTri,:) ;
        nElems = size(Elems,1) ;
        if nElems==0 ; return ; end
        
    % Points: create "satellite" points defining a polygon
        isPoint = sum(Elems>0,2)==1 ;
        if any(isPoint)
        % Get involved node indices and reshape the element list
            nn = Elems(isPoint,1) ;
            Elems = Elems(~isPoint,:) ; % remove point elements
        % Define the satellite points
            nDup = this.PointSatellites ; % if==4, create a square
            theta = 2*pi/nDup*((0:nDup-1)'+.5) ; % satellite angles
            S = repmat(Vertices(nn,:),[nDup 1]) + repelem(this.EdgeThickness/2*[cos(theta) sin(theta)],numel(nn),1) ;
            Vertices = [Vertices ; S] ;
        % Transfer matrix
            N2Q = sparse(1:numel(nn),nn(:)',1,numel(nn),nNodes) ;
            N2Q = repmat(N2Q,[nDup 1]) ;
            N2V = kron(speye(2),[speye(nNodes);N2Q])*N2V ;
        % New Elements
            if size(Elems,2)<nDup ; Elems(:,end:nDup) = NaN ; end
            newElems = nNodes + (1:numel(nn))' + numel(nn)*(0:nDup-1) ;
            Elems = [Elems ; newElems] ;
        end
        
    % Lines: create a quadrilateron from offsetted lines
        isLine = sum(Elems>0,2)==2 ;
        if any(isLine)
            mergeAdj = true ; % merge adjacent nodes for strain regularisation ?
            nodeNormals = false ; % offset with node normals ?
        % Get involved node indices and reshape the element list
            nn = Elems(isLine,1:2) ;
            Elems = Elems(~isLine,:) ; % remove line elements
            nLines = size(nn,1) ; 
            if mergeAdj 
                [un,inn,iun] = unique(nn(:)) ; 
                Mu = sparse(iun(:)',1:numel(nn),1,numel(un),numel(nn)) ;
                Mu = diag(1./sum(Mu,2))*Mu ; % mean over duplicated nodes
            end
        % Line to Nodes transfer for offsetting
            nNodes = size(Vertices,1) ;
            if nodeNormals
                l2n = sparse(nn(:)',repmat(1:nLines,[1 2]),1,nNodes,nLines) ;
                l2n = diag(1./sum(l2n,2))*l2n ; % mean over adjacent edges
                l2n = l2n(nn(:),:) ;
            else
                l2n = sparse(1:2*nLines,repmat(1:nLines,[1 2]),1,2*nLines,nLines) ;
            end
        % Define the offset points
            h = this.EdgeThickness/2 ; % position shift along normal
            RotMat = sparse([0 1 ; -1 0]) ; % 90°rotation matrix tangent->normal
            tan = Vertices(nn(:,2),:)-Vertices(nn(:,1),:) ; % line tangent
            nor = tan*RotMat ; % line normal
            Pnor = l2n*nor ; % normal on line nodes for offset
            P = h*Pnor./(sqrt(sum(Pnor.^2,2))) ; % offset distance
            P = repmat(Vertices(nn,:),[2 1]) + [-P ; P] ;
            if mergeAdj ; P = blkdiag(Mu,Mu)*P ; end
            Vertices = [Vertices ; P] ;
        % Transfer matrix
            % Translation part: dp_t = u
                N2Qt = sparse(1:numel(nn),nn(:)',1,numel(nn),nNodes) ;
                N2Qt = repmat(N2Qt,[2 1]) ;
                N2Qt = kron(speye(2),N2Qt) ;
            % Rotation part: dp_r = h/L*(RotMat-kron(nor,tan))*D(u)
                L = sqrt(sum(tan.^2,2)) ; 
                e = diag(sparse(1./L))*h ; % excentricity
                E = [e;e;-e;-e] ; % excentricities of all new vertices
                D = sparse(1:nLines,nn(:,2),1,nLines,nNodes)-sparse(1:nLines,nn(:,1),1,nLines,nNodes) ; % difference matrix
                N2Qr = kron(RotMat,E*D) ; % h/L*RotMat*D(u)
                nt = (L.^-2).*nor.*permute(tan,[1 3 2]) ; % kron(nor,tan) [nLines 2 2]
                ent = sparse(E*nt(:,:)) ; % h/L*kron(nor,tan) [4*nLines 4] ;
                N2Qr = N2Qr + [diag(ent(:,1)) diag(ent(:,3)) ; diag(ent(:,2)) diag(ent(:,4))]*kron(speye(2),repmat(D,[4 1])) ;
            % Total matrix recursion
                N2Q = N2Qt+N2Qr ; 
                if mergeAdj ; N2Q = blkdiag(Mu,Mu,Mu,Mu)*N2Q ; end
                I = speye(nNodes) ; O = sparse(nNodes,nNodes) ;
                N2Q = [I O ; N2Q(1:end/2,:) ; O I ; N2Q(end/2+1:end,:)] ;
                N2V = N2Q*N2V ;
        % New Elements
            if size(Elems,2)<4 ; Elems(:,end:4) = NaN ; end
            newElems = (1:nLines)' + nLines*([0 1 3 2]) ;
            if mergeAdj 
                newInd = [iun;iun+numel(un)] ; 
                newElems = newInd(newElems) ; 
            end
            newElems = newElems + nNodes ;
            Elems = [Elems ; newElems] ;
        end
        
        
    % Polygonal elements: divide by simplices defined by edges and centroid
        nNodes = size(Vertices,1) ;
        nElems = size(Elems,1) ;
        % Transfer from nodes to elem centroids
            ee = repmat((1:nElems)',[1 size(Elems,2)]) ;
            nn = Elems ;
            valid = ~isnan(Elems) ;
            N2C = sparse(ee(valid),nn(valid),1,nElems,nNodes) ;
            N2C = diag(1./sum(N2C,2))*N2C ; % mean over element nodes
        % Add centroids to vertices
            N2V = kron(speye(2),[speye(nNodes);N2C])*N2V ; 
            Vertices = [Vertices ; N2C*Vertices] ;
        % Edge node indices
            e1 = Elems ;
            Nn = sum(~isnan(e1),2) ;
            e2 = circshift(e1,-1,2) ;
            e2(sub2ind(size(e2),1:size(e2,1),Nn')) = e2(:,end) ;
        % Triangles
            newTris = [e1(:) e2(:) repmat(nNodes+(1:nElems)',[size(Elems,2) 1])] ;
            newTris = newTris(all(newTris>0,2),:) ; % keep valid simplices
        % Add to the triangle list
            Tris = [Tris ; newTris] ;
            
    % Delete points that are not associated to any simplex
        lonePt = ~ismember(1:size(Vertices,1),Tris(:)) ;
        Vertices(lonePt,:) = [] ;
        N2V([lonePt lonePt],:) = [] ;
        newInd = [NaN cumsum(~lonePt)] ;
        Tris(isnan(Tris)) = 0 ;
        Tris = newInd(Tris+1) ;
        
    end
    
    function pl = plotMesh(this)
    % Plot the underlying simplex mesh for debugging purposes
        axis equal tight ij
        imagesc(this.RefImgs{end}) ;
        colormap gray
        m = pkg.geometry.mesh.Mesh('Nodes',this.Vertices,'Elems',this.Tris) ;
        pl = plot(m,'FaceColor','none','EdgeColor','r') ;
    end
    
    function [MAP,IN] = P1ShapeFunctions(this,refImg,Vertices,Tris)
    % Compute the value of P1 triangular element shape functions at pixel
    % coordinates of the reference image
        if nargin<2 ; refImg = this.RefImgs{end}(:,:,1) ; end % gray-scale reference
        if nargin<3 ; Vertices = this.Vertices ; end % reference nodes
        if nargin<4 ; Tris = this.Tris ; end % seed elements
        
    % Infos
        nVertices = size(Vertices,1) ;
        nTris = size(Tris,1) ;
        [nI,nJ] = size(refImg,[1 2]) ;

    % 1) In which element bounding box is which pixel ?
    % Element bounding boxes
        xe = Vertices(Tris(:),:) ; % [nTris*3 2]
        xe = reshape(xe,[nTris 3 2]) ; % [nTris 3 2] 
        xmin = min(floor(xe),[],2) ; % [nTris 1 2]
        xmax = max(ceil(xe),[],2) ; % [nTris 1 2]
    % concatenate all pixels
        nPixInBBoxDim = xmax-xmin+1 ; % [nTris 1 2]
        nPixInBBox = prod(nPixInBBoxDim,3) ; % [nTris 1 1]
        nPixToTest = sum(nPixInBBox,1) ; % [1 1 1]
        ee = repelem(1:nTris,nPixInBBox) ; % element index
        nPixInPreviousBBoxes = [0 cumsum(nPixInBBox(1:end-1)')] ;
        pp = (0:nPixToTest-1) - repelem(nPixInPreviousBBoxes,nPixInBBox) ; % pixel index
        ii = mod(pp,repelem(nPixInBBoxDim(:,2)',nPixInBBox)) ; % bbox row index
        jj = (pp-ii)./repelem(nPixInBBoxDim(:,2)',nPixInBBox) ; % bbox column index
        ii = ii + xmin(ee,2)' ; % image row index
        jj = jj + xmin(ee,1)' ; % image column index

    % 2) Compute the shape functions
    % Relative pixel coordinates
        PP = [jj(:) ii(:)] ; % Pixel coordinates [nPixToTest 2]
        xe = permute(xe,[1 3 2]) ; % triangle nodes [nTris 2 3]
        PP = PP-xe(ee,:,1) ; % Relative pixel coordinates [nPixToTest 2]
    % Local coordinates xi such that A.xi = PP so xi = inv(A).PP
        A = xe(:,:,[2 3])-xe(:,:,1) ; % edge vectors [nTris 2 2]
        detA = A(:,1,1).*A(:,2,2)-A(:,2,1).*A(:,1,2) ; % [nTris 1 1]
        invA = [A(:,2,2) -A(:,2,1) -A(:,1,2) A(:,1,1)]./detA ; % [nTris 4]
        invA = reshape(invA,[nTris 2 2]) ; % [nTris 2 2] ;
        xi = sum(invA(ee,:,:).*reshape(PP,[nPixToTest 1 2]),3) ; % [nPixToTest 2]
    % Shape Functions N = [1-Xi1-Xi2 Xi1 Xi2]
        NN = [1-xi(:,1)-xi(:,2) xi(:,1) xi(:,2)] ;
        
    % 3) row-col to pixel indices
        nPix = nI*nJ ;
        pp = ii+(jj-1)*nI ;

    % 4) keep only valid tests (so that shape functions are in [0 1])
        valid = all(NN>=0 & NN<=1,2) ;
        pp = pp(valid) ;
        NN = NN(valid,:) ;
        ee = ee(valid) ;
        
    % 5) attach only one element by pixel (can happen at edges)
       [pp,iu] = unique(pp) ;
       NN = NN(iu,:) ;
       ee = ee(iu) ;

    % 6) Build the sparse matrix
    % Node indices
        nn = Tris(ee,:) ;
    % Sparse
        MAP = sparse(repmat(pp(:),[3 1]),nn(:),NN(:),nPix,nVertices) ;
        
    % 4) build the INSIDE matrix
        if nargout>=2 ; IN = sparse(pp,ee,true,nPix,nTris) ; end
    end
    
    function [dB,Hr] = strainGradientRegularisation(this,normalize)
    % Strain gradient regularization matrices
        if nargin<2 ; normalize = true ; end
    % Create a temporary seed with the augmented simplex mesh
        seed = copy(this.Seed) ; seed.Elems = this.Tris ; seed.Points = this.Vertices ;
        Elems = seed.Elems ; nElems = size(Elems,1) ;
        Nodes = seed.Points ; nNodes = size(Nodes,1) ;
    % Gradient matrices D1 & D2
        [D1,D2] = seed.diffMat(Nodes,false) ;
    % Strain matrix B so that [E11(:) ; E22(:) ; 2*E12(:)] = B*[U]
        O = sparse(nElems,nNodes) ; % matrix full of zeros
        B = [D1 O ; O D2 ; D2 D1] ;
    % Edge differentiation matrix G so that df = G*[f] where f is defined on
    % elements and df measures the difference between 2 element values 
        % element->edge connectivity
        [~,ele2edg] = seed.getEdges() ;
        % Keep interior edges only 
        boundaryEdges = sum(ele2edg,2)==1 ;
        ele2edg(boundaryEdges,:) = [] ;
        % Difference matrix
        [elem,edg] = find(ele2edg') ;
        val = repmat([1;-1],[numel(edg)/2 1]) ;
        G = sparse(edg,elem,val,size(ele2edg,1),nElems) ;
    % Distance between element centroids
        xe = reshape(Nodes(Elems(:),:),[nElems 3 2]) ;
        C = reshape(mean(xe,2),[nElems 2]) ; % element centroids
        dC = sqrt(sum((G*C).^2,2)) ; % distance between centroids
    % Strain Gap
        if normalize ; G = (1./dC(:)).*G ; end % difference normalized by the centroid distance
        dB = blkdiag(G,G,G)*B ; % dB*[U] measures the (normalized) strain gap
    % Weighting with mean element area (surface integral)
        w = polyarea(xe(:,:,1),xe(:,:,2),2) ;
        w = w(:)./max(sum(ele2edg,1),eps)' ;
        W = diag(double(logical(G))*sparse(w)) ;
        W = blkdiag(W,W,W) ;
    % Associated Hessian
        if nargout>=2 ; Hr = dB'*W*dB ; Hr = (Hr+Hr')/2 ; end
        dB = W*dB ;
    end
    
    function val = bilinearInterp(this,img,xx,extrapVal)
    % Bilinear interpolation of an image at real-valued coordinates xx
    % faster than gg = interp2(img,xx(:,1),xx(:,2),'linear') ;
    % xx = [nValues 2] = [jj(:) ii(:)] ;
    % with ii and jj resp. the (real-valued) row and column indices
    % extrapVal: extrapolation value (default: NaN)

        if nargin<4 ; extrapVal = NaN ; end

        % Image informations
        [nI,nJ,~] = size(img) ;

        % Valid coordinates
        valid = xx(:,1)<=nJ ...
                & xx(:,1)>=1 ...
                & xx(:,2)<=nI ...
                & xx(:,2)>=1 ;

        % Dummy values
        xx(~valid,:) = 1 ;

        % Integer part of the cordinates
        ji = floor(xx) ;
        ji(ji(:,1)==nJ,1) = nJ-1 ;
        ji(ji(:,2)==nI,2) = nI-1 ;

        % Neightboring pixels
        p1 = ji(:,2)+nI*(ji(:,1)-1) ;
        p2 = p1 + nI ; 
        p3 = p2+1 ; 
        p4 = p1+1 ;

        % residual coordinates
        dx = xx-ji ;

        % bilinear interpolation
        val = img(p1).*(1-dx(:,1)).*(1-dx(:,2)) ...
            + img(p2).*dx(:,1).*(1-dx(:,2)) ...
            + img(p3).*dx(:,1).*dx(:,2) ...
            + img(p4).*(1-dx(:,1)).*dx(:,2) ;

        % Non-valid values
        val(~valid) = extrapVal ;

    end
end


end

