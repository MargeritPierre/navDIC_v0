classdef navDIC_GlobalDICMacro < navDIC_DICMacro
%NAVDIC_GlobalDICMacro Implement Global DIC

properties
    % Image interpolation
        InterpOrder = 'linear' ;
    % Image difference criterion
        DiffCriterion = ... 'diff' ... Simple difference
                        ... 'zeromean' ... Zero-mean difference
                         'normalized' ... Normalized Zero-mean difference
                         ;
    % Descent Algorithm
        GNAlgorithm = ... 'full' ... full Gauss-Newton
                  'modified' ... modified (assume "grad(g(x+u))=grad(f(x))") OK for small perturbations
                    ;
    % Geometry enhancement
        EdgeThickness = 30 ; % include surrounding pixels of (outer) edges
        PointSatellites = 4 ; % add satellite points to isolated points
    % Geometry validation criteria
        CullOutOfFrame = true ; % Cull out of frame points
        MinCorrCoeff = 0 ; % Below this, elements are culled
        MaxMeanElemResidue = Inf ; % Above this, elements are culled
        ThresholdValidGeometry = 0 ; % Check correlation. coeffs when the (normA/minNorm)<thresholdValidGeometry. (0 disable the check)
    % Regularisation
        StepRatio = 1 ; % Descent step ratio, damping the convergence is <1
        RegCrit = 'rel' ; % second gradient minimization: absolute variation ('abs') or relative ('rel')
        Beta = 1*1e7 ; % Strain gradient penalisation coefficient
        EpsTrsh = 1e0 ; % Limit value for the regularisation weights (active when regCrit = 'rel')
    % Convergence Criteria
        MaxIt = 100 ; % Maximum number of Newton-Raphson iterations
        MaxDisp = 1e-4 ; % Maximum displacement of a node
        MaxResidueRelativeVariation = Inf ; % Maximum relative variation of the image residue (RMSE)
        MinCorrdU = -.999 ; % Maximum update correlation between two consecutive iterations (avoids oscillations)
    % Plotting
        Debug logical = false
        PlotRate = Inf ; % Plot Refresh Frequency 
        PlotEachIteration = false ; % Plot at every iteration (without necessary pausing, bypass plotRate)
        PlotEachFrame = false ; % Plot at every Frame end (without necessary pausing, bypass plotRate)
        PauseAtPlot = false ; % Pause at each iteration for debugging
end

methods
    function this = navDIC_GlobalDICMacro()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setupUI(this,hd)
    % Setup the macro (change parameters, etc)
        hd = setupUI@navDIC_DICMacro(this,hd) ; % Set parameters common to DIC macros
        if isempty(this.Seed) ; return ; end
    % Set globalDIC parameters
        defInputs = { ...
                        'Gauss-Newton descent algorithm [full/modified]' , this.GNAlgorithm ...
                        ; 'Image difference criterion [diff/zeromean/normalized]' , this.DiffCriterion ...
                        ; 'Maximum displacement update at convergence' , num2str(this.MaxDisp) ...
                        ; 'Strain regularisation parameter Beta' , num2str(this.Beta) ...
                        ; 'Strain regularisation mode [abs/rel]' , this.RegCrit ...
                        ; 'Maximum number of iterations' , num2str(this.MaxIt) ...
                        ; 'Step ratio [0->1]' , num2str(this.StepRatio) ...
                        ; 'Image interpolation order [linear/cubic]' , num2str(this.InterpOrder) ...
                        ; 'Edge Thickness' , num2str(this.EdgeThickness) ...
                        ; 'Number of satellites for isolated points' , num2str(this.PointSatellites) ...
                        ; 'Debug mode' , num2str(this.Debug) ...
                    } ;
        out = inputdlg(defInputs(:,1),'Global DIC parameters',1,defInputs(:,2)) ;
        if isempty(out) ; return ; end
        this.GNAlgorithm = out{1} ;
        this.DiffCriterion = out{2} ; 
        this.MaxDisp = str2double(out{3}) ;
        this.Beta = str2double(out{4}) ;
        this.RegCrit = out{5} ;
        this.MaxIt = str2double(out{6}) ;
        this.StepRatio = str2double(out{7}) ;
        this.InterpOrder = str2double(out{8}) ;
        this.EdgeThickness = str2double(out{9}) ;
        this.PointSatellites = str2double(out{10}) ;
        this.Debug = str2double(out{11}) ;
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
properties
   Vertices % (Augmented) simplex mesh nodes [nVertices 2]
   Tris % (Augmented) simplicex mesh elements [nTris 3]
   N2V % transfer matrix from the seed mesh nodes to the simplex mesh vertices [2*nVertices 2*nNodes]
   MAP % Sparse nodal->pixels interpolation matrix [nPixelsInROI nNodes] (FEM shape functions evaluated at pixel coordinates)
   IN % Sparse logical matrix [nPixelsInROI nElems]: which pixel is in which element ?
   ROI % Region Of Interest (everywhere the mapping is not zero)
   Grad % Sparse matrices aiming to compute the image gradient {[nPixelsInROI nPixels]}(nCoord)
   dr_da % derivative of the image difference w.r.t parameters
   H % hessian
   Hr % regularization hessian
end
methods
    function setupDIC(this,hd)
    % Prepare the DIC data
        refImg = this.RefImgs{end}(:,:,1) ; % gray-scale reference
        [nI,nJ] = size(refImg,[1 2]) ;
    % Build the augmented simplex mesh
        [this.Vertices,this.Tris,this.N2V] = this.simplexMesh() ;
    % Construct the FEM interpolation mapping
        [this.MAP,this.IN] = this.P1ShapeFunctions(refImg) ;
    % Get the ROI
        this.ROI = reshape(any(this.IN,2),[nI nJ]) ;
        this.ROI([1 end],:) = false ; % remove image borders
        this.ROI(:,[1 end]) = false ; % remove image borders
        this.MAP = this.MAP(this.ROI,:) ; % reduce the problem to ROI
        this.IN = this.IN(this.ROI,:) ; % reduce the problem to ROI
    % Compute regularisation data
        [~,this.Hr] = this.strainGradientRegularisation() ;
        this.Hr = this.N2V'*this.Hr ;
    % Compute the GN data
        [this.dr_da,this.H] = this.GNData(refImg) ;
    end
    
    function X = updateDIC(this,X,imgs,refImgs)
    % Update the configuration X using DIC performed on imgs
    % For now, only working on one gray-scale image
        G = refImgs{end}(:,:,1) ; % reference <MONOCHROME SINGLE CAMERA!>
        g = imgs{end}(:,:,1) ; % current <MONOCHROME SINGLE CAMERA!>
    % Prepare the reference frame
        G = G(this.ROI) ;
        G = this.normalizedImg(G) ;
    % Initialize
        U = X*0 ; it = 0 ; outFlag = [] ;
        dx = this.Vertices-reshape(this.N2V*this.Seed.Points(:),[],2) ;
    % DIC Loop
        while isempty(outFlag)
        % Image interpolation
            x = reshape(this.N2V*(X(:)+U(:)),[],2) ; 
            xx = this.MAP*(x+dx) ; % [nPix nCoord]
            gx = this.bilinearInterp(g,xx) ; % [nPix 1]
        % Image residual
            r = this.normalizedImg(gx(:))-G(:) ;
        % Full gauss-newton ?
            switch this.GNAlgorithm
                case 'modified'
                otherwise
                    [this.dr_da,this.H] = this.GNData(g,xx) ;
            end
        % Jacobians
            j = this.dr_da'*r ;
            switch this.RegCrit
                case 'abs' ; jr = this.Hr*(dx(:)+this.N2V*(X(:)+U(:))) ;
                case 'rel' ; jr = this.Hr*(this.N2V*U(:)) ;
            end % full H & j
            j = j + this.Beta*jr ;
            H = this.H + this.Beta*this.Hr*this.N2V ;
        % Select components
            switch this.DispComp
                case 'both' % do nothing
                case 'X'
                    j = j(1:end/2) ;
                    H = H(1:end/2,1:end/2) ;
                case 'Y'
                    j = j(end/2+1:end) ;
                    H = H(end/2+1:end,end/2+1:end) ;
            end
        % Update
            du = - H \ j ; 
            switch this.DispComp
                case 'both' % do nothing
                    dU = reshape(du,[],2) ;
                case 'X'
                    dU = [1 0].*du ;
                case 'Y'
                    dU = [0 1].*du ;
            end
            U = U + dU ;
        % Break criterion
            it = it+1 ;
            if it>=this.MaxIt ; outFlag = 'maximum iterations reached' ; end
            if max(abs(dU))<=this.MaxDisp ; outFlag = 'below displacement threshold' ; end
            if all(isnan(dU)) ; outFlag = 'all nodes are NaN' ; end
        % Display infos
            disp([ ...
                this.Name ':' ...
                ' it ' num2str(it) ...
                ' max(dU) ' num2str(max(abs(dU))) ...
                ]) ;
        % Debug ?
            if this.Debug
                global hd
            % Update seed data
                this.Seed.MovingPoints(:,:,hd.CurrentFrame) = X+U ;
                this.Seed.computeDataFields([],hd.CurrentFrame) ;
            % Update navDIC previews
                hd = updateAllPreviews(hd) ;
            % Draw
                drawnow ;
            end
        end
    % Return the config
        X = X+U ;
    end
    
    function [dr_da,H] = GNData(this,img,xx)
    % Return Gauss-Newton residual derivative and Hessian
    % Compute the image gradient
        dI_dx = conv2(img,[1 , 0 , -1]/2,'same') ;
        dI_dy = conv2(img,[1 ; 0 ; -1]/2,'same') ;
    % Interpolate at the given coordinates & compute the residual
        if nargin<3 % take the ROI by default
        %dr_da = [diag(sparse(this.Grad{1}*img(:)))*this.MAP diag(sparse(this.Grad{2}*img(:)))*this.MAP] ;
            imgdata = img(this.ROI) ;
            dr_da = [diag(sparse(dI_dx(this.ROI)))*this.MAP diag(sparse(dI_dy(this.ROI)))*this.MAP] ;
        else % 
            imgdata = this.bilinearInterp(img,xx) ;
            dI_dx = this.bilinearInterp(dI_dx,xx) ;
            dI_dy = this.bilinearInterp(dI_dy,xx) ;
            dr_da = [diag(sparse(dI_dx))*this.MAP diag(sparse(dI_dy))*this.MAP] ;
        end
    % Project to seed nodes
        dr_da = dr_da*this.N2V ;
    % Normalization weights
        [~,W] = this.normalizedImg(imgdata) ;
    % Hessian
        H = dr_da'*diag(sparse(W(:)))*dr_da ;
        H = (H+H')/2 ;
    end
    
    function [I,W] = normalizedImg(this,I,w)
    % Return the normalized image
        W = ones(size(I)) ; % weighting
        if ismember(this.DiffCriterion,{'diff'}) ; return ; end
    % Weights
        if nargin<3 ; w = this.IN ; end
        sumW = sum(w,1)' ;
    % Substract the mean image
        meanI = (w'*I)./sumW ;
        I = I-w*meanI ;
        if ismember(this.DiffCriterion,{'zeromean'}) ; return ; end
    % Normalize
        normI = sqrt( ( w'*(I.^2) )./sumW ) ;
        W = w*(1./normI) ;
        I = I.*W ;
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
end


end

