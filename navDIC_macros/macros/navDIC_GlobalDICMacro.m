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
    % Geometry validation criteria
        CullOutOfFrame = true ; % Cull out of frame points
        MinCorrCoeff = 0 ; % Below this, elements are culled
        MaxMeanElemResidue = Inf ; % Above this, elements are culled
        ThresholdValidGeometry = 0 ; % Check correlation. coeffs when the (normA/minNorm)<thresholdValidGeometry. (0 disable the check)
    % Regularisation
        StepRatio = 1 ; % Descent step ratio, damping the convergence is <1
        RegCrit = 'abs' ; % second gradient minimization: absolute variation ('abs') or relative ('rel')
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
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
        hd = setup@navDIC_DICMacro(this,hd) ; % Set parameters common to DIC macros
        if isempty(this.Seed) ; return ; end
    % Set globalDIC parameters
        defInputs = { ...
                        'Gauss-Newton descent algorithm [full/modified]' , this.GNAlgorithm ...
                        ; 'Image difference criterion [diff/zeromean/normalized]' , this.DiffCriterion ...
                        ; 'Maximum displacement update at convergence' , num2str(this.MaxDisp) ...
                        ; 'Strain regularisation parameter Beta' , num2str(this.Beta) ...
                        ; 'Maximum number of iterations' , num2str(this.MaxIt) ...
                        ; 'Step ratio [0->1]' , num2str(this.StepRatio) ...
                        ; 'Image interpolation order [linear/cubic]' , num2str(this.InterpOrder) ...
                        ; 'Debug mode' , num2str(this.Debug) ...
                    } ;
        out = inputdlg(defInputs(:,1),'Global DIC parameters',1,defInputs(:,2)) ;
        if isempty(out) ; return ; end
        this.GNAlgorithm = out{1} ;
        this.DiffCriterion = out{2} ; 
        this.MaxDisp = str2double(out{3}) ;
        this.Beta = str2double(out{4}) ;
        this.MaxIt = str2double(out{5}) ;
        this.StepRatio = str2double(out{6}) ;
        this.InterpOrder = str2double(out{7}) ;
        this.Debug = str2double(out{8}) ;
    % Prepare the DIC data
        this.setupDIC(hd) ;
    end
end


%% DIC PROCEDURES
properties
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
    % Construct the FEM interpolation mapping
        [this.MAP,this.IN] = this.P1ShapeFunctions(refImg) ;
    % Get the ROI
        this.ROI = reshape(any(this.IN,2),[nI nJ]) ;
        this.MAP = this.MAP(this.ROI,:) ;
        this.IN = this.IN(this.ROI,:) ;
    % Construct the gradient matrices
        pp = find(this.ROI(:)) ;
        ip = (1:numel(pp))' ;
        d_dx1 = sparse(ip.*[1 1],pp+[-1 1]*nI,pp*0+[-1 1]/2,numel(pp),nI*nJ) ;
        d_dx2 = sparse(ip.*[1 1],pp+[-1 1],pp*0+[-1 1]/2,numel(pp),nI*nJ) ;
        this.Grad = {d_dx1,d_dx2} ;
    % Compute the GN data
        [this.dr_da,this.H] = GNData(this,refImg) ;
    % Compute regularisation data
        [~,this.Hr] = strainGradientRegularisation(this) ;
    end
    
    function X = updateDIC(this,X,imgs)
    % Update the configuration X using DIC performed on imgs
    % For now, only working on one gray-scale image
        G = this.RefImgs{end}(:,:,1) ; % reference
        g = imgs{end}(:,:,1) ; % current
    % Prepare the reference frame
        G = G(this.ROI) ;
        G = this.normalizedImg(G) ;
    % Initialize
        U = X*0 ; it = 0 ; outFlag = [] ;
    % DIC Loop
        while isempty(outFlag)
        % Image interpolation
            xx = this.MAP*(X+U) ; % [nPix nCoord]
            gx = this.bilinearInterp(g,xx) ; % [nPix 1]
        % Image residual
            r = this.normalizedImg(gx(:))-G(:) ;
        % Jacobians
            j = this.dr_da'*r ;
            switch this.RegCrit
                case 'abs' ; jr = this.Hr*(X(:)+U(:)) ;
                case 'rel' ; jr = this.Hr*U(:) ;
            end
        % Update
            dU = - (this.H + this.Beta*this.Hr) \ (j + this.Beta*jr) ;
            U = U + reshape(dU,[],2) ;
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
                this.Seed.MovingPoints(:,:,hd.CurrentFrame) = X ;
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
    
    function [dr_da,H] = GNData(this,img)
    % Return Gauss-Newton residual derivative and Hessian
        dr_da = [(this.Grad{1}*img(:)).*this.MAP (this.Grad{2}*img(:)).*this.MAP] ;
        [~,W] = this.normalizedImg(img(this.ROI)) ;
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
    
    function I = transformImage(this,X,x,I)
    % Transform the image(s) I from configuration X to x
        I(this.ROI) = this.bilinearInterp(I,this.MAP*x) ;
    end
end

%% UTIL FUNCTIONS
methods
    function [MAP,IN] = P1ShapeFunctions(this,refImg,Nodes,Elems)
    % Compute the value of P1 triangular element shape functions at pixel
    % coordinates of the reference image
        if nargin<2 ; refImg = this.RefImgs{end}(:,:,1) ; end % gray-scale reference
        if nargin<3 ; Nodes = this.Seed.Points ; end % reference nodes
        if nargin<4 ; Elems = this.Seed.Elems(:,1:3) ; end % triangular elements
    % Infos
        nNodes = size(Nodes,1) ;
        nElems = size(Elems,1) ;
        [nI,nJ] = size(refImg,[1 2]) ;

    % 1) In which element bounding box is which pixel ?
    % Element bounding boxes
        xe = Nodes(Elems(:),:) ; % [nElems*3 2]
        xe = reshape(xe,[nElems 3 2]) ; % [nElems 3 2] 
        xmin = min(floor(xe),[],2) ; % [nElems 1 2]
        xmax = max(ceil(xe),[],2) ; % [nElems 1 2]
    % concatenate all pixels
        nPixInBBoxDim = xmax-xmin+1 ; % [nElems 1 2]
        nPixInBBox = prod(nPixInBBoxDim,3) ; % [nElems 1 1]
        nPixToTest = sum(nPixInBBox,1) ; % [1 1 1]
        ee = repelem(1:nElems,nPixInBBox) ; % element index
        nPixInPreviousBBoxes = [0 cumsum(nPixInBBox(1:end-1)')] ;
        pp = (0:nPixToTest-1) - repelem(nPixInPreviousBBoxes,nPixInBBox) ; % pixel index
        ii = mod(pp,repelem(nPixInBBoxDim(:,2)',nPixInBBox)) ; % bbox row index
        jj = (pp-ii)./repelem(nPixInBBoxDim(:,2)',nPixInBBox) ; % bbox column index
        ii = ii + xmin(ee,2)' ; % image row index
        jj = jj + xmin(ee,1)' ; % image column index

    % 2) Compute the shape functions
    % Relative pixel coordinates
        PP = [jj(:) ii(:)] ; % Pixel coordinates [nPixToTest 2]
        xe = permute(xe,[1 3 2]) ; % triangle nodes [nElems 2 3]
        PP = PP-xe(ee,:,1) ; % Relative pixel coordinates [nPixToTest 2]
    % Local coordinates xi such that A.xi = PP so xi = inv(A).PP
        A = xe(:,:,[2 3])-xe(:,:,1) ; % edge vectors [nElems 2 2]
        detA = A(:,1,1).*A(:,2,2)-A(:,2,1).*A(:,1,2) ; % [nElems 1 1]
        invA = [A(:,2,2) -A(:,2,1) -A(:,1,2) A(:,1,1)]./detA ; % [nElems 4]
        invA = reshape(invA,[nElems 2 2]) ; % [nElems 2 2] ;
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
        nn = Elems(ee,:) ;
    % Sparse
        MAP = sparse(repmat(pp(:),[3 1]),nn(:),NN(:),nPix,nNodes) ;
        
    % 4) build the INSIDE matrix
        if nargout>=2 ; IN = sparse(pp,ee,true,nPix,nElems) ; end
    end
    
    function [dB,Hr] = strainGradientRegularisation(this,normalize)
    % Strain gradient regularization matrices
        if nargin<2 ; normalize = true ; end
        Elems = this.Seed.Elems ; nElems = size(Elems,1) ;
        Nodes = this.Seed.Points ; nNodes = size(Nodes,1) ;
    % Gradient matrices D1 & D2
        [D1,D2] = this.Seed.diffMat(Nodes,false) ;
    % Strain matrix B so that [E11(:) ; E22(:) ; 2*E12(:)] = B*[U]
        O = sparse(nElems,nNodes) ; % matrix full of zeros
        B = [D1 O ; O D2 ; D2 D1] ;
    % Edge differentiation matrix G so that df = G*[f] where f is defined on
    % elements and df measures the difference between 2 element values 
        % element->edge connectivity
        [~,ele2edg] = this.Seed.getEdges() ;
        % Keep interior edges only 
        boundaryEdges = sum(ele2edg,2)==1 ;
        ele2edg(boundaryEdges,:) = [] ;
        % Difference matrix
        [elem,edg] = find(ele2edg') ;
        val = repmat([1;-1],[numel(edg)/2 1]) ;
        G = sparse(edg,elem,val,size(ele2edg,1),nElems) ;
    % Distance between element centroids
        xe = reshape(Nodes(Elems(:),:),[nElems 3 2]) ;
        C = reshape(mean(xe,2),[nElems 2]) ;
        dC = sqrt(sum((G*C).^2,2)) ;
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

