%% GLOBAL DIC USING THE PKG.GEOMETRY.MESH
% Generalizes the current GlobalDIC implementation

%% INITIALIZATION
% CLEAN WORKSPACE
    clc
    global hd
    clearvars -except hd

% INITIALIZATION PARAMETERS
    camID = 1 ;
    seedNumber = 1 ;
    frames = '[1:end]' ; % Frames taken for DIC (allows decimation)
    dicDir = -1 ; % DIC running direction ('forward=1' or 'backward=-1')
    refFrame = 'last' ; % Reference image ('first' , 'last' or number)
    refConfig = 'Nodes' ; % Reference configuration: 'Nodes' (as meshed) or 'Current' (uses preceding computed displacement)
    strainCriterion = 'full' ;
    averagePreviousFrames = true ; % Ref frame is the average of the previous/next ones in forward/backward modes
    normToImageClassRange = true ; % Normalize images to their dataclass range
    timeMeanLength = 0 ; % Time averaging of images
    showInit = true ;
    codeProfile = false ; % Code timing
    figTag = 'Global DIC' ;

% Load Frames
    globalDIC_01_LoadFrames ;
    
% Process the Seed
    seed = hd.Seeds(seedNumber) ;
    Nodes = seed.Points ; nNodes = size(Nodes,1) ;
    Elems = seed.Elems ; nElems = size(Elems,1) ;
    dicMesh = pkg.geometry.mesh.Mesh('Nodes',Nodes,'Elems',Elems) ;
    Edges = dicMesh.Edges.NodeIdx ; nEdges = dicMesh.nEdges ;
    
% Mesh connectivity
    tri2nod = dicMesh.elem2node ;
    nNodesByElements = dicMesh.Elems.nNodes ;
    edg2nod = dicMesh.edge2node ;
    tri2edg = dicMesh.elem2edge ;
    
% Display Initialization
    if showInit ; globalDIC_03_1_ShowInitialization ; end

    
% SHAPE FUNCTIONS

% Parameters
    pixTol = 1e-4 ; % pixel localization tolerance

% Pixel localization (in local element coordinates)
% Which pixel is in which element bounding box ?
    xe = dicMesh.Elems.dataAtIndices(dicMesh.Nodes) ;
    bbox = [min(xe,[],2) max(xe,[],2)] ;
    [ji,ee] = pkg.data.domainIndices(bbox(:,:,:)) ;
% Localize
    [E,ie] = dicMesh.localize(ji,dicMesh.Elems... localize in elements
                            ,false ... without extrapolation
                            ,dicMesh.Nodes ... in the reference config
                            ,pixTol ... with reduced tolerance
                            ,(1:numel(ee))' ... for all points
                            ,ee ... with respect to previously found elements
                            ) ;
% Cull invalid elements
    valid = all(~isnan(E),2) ;
    E = E(valid,:) ;
    ie = ie(valid) ;
    ROI = ji(valid,:) ;
    nROI = size(ROI,1) ;
% Sort ROI (for memory access)
    ind = max(ROI(:,2))*(ROI(:,1)-1) + ROI(:,1) ; 
    [~,is] = sort(ind) ;
    ROI = ROI(is,:) ;
    E = E(is,:) ;
    ie = ie(is) ;
% Correlation domain
    indDOMAIN = ROI(:,2) + nI*(ROI(:,1)-1) ;
    IId = II(indDOMAIN) ;
    JJd = JJ(indDOMAIN) ;
% Sparse matrix of shape functions
    MAPPING = dicMesh.interpMat(E,ie) ; % [nROI nNodes]
% Sparse (boolean) matrix of pixels in elements
    INSIDE = sparse(1:nROI,ie(:)',true,nROI,dicMesh.nElems) ;
    
    
% SECOND GRADIENT CRITERION: 

% EDGES (all elements)

    % Gradient matrix
        Gr = dicMesh.gradMat ;
        O = sparse(dicMesh.nElems,dicMesh.nNodes) ;
        G = [Gr{1} O ; O Gr{1} ; Gr{2} O ; O Gr{2}] ;

    % Linearized Green-Lagrange strains
        B = [Gr{1} O ; O Gr{2} ; Gr{2} Gr{1}] ;

    % Edge differenciation matrix
            bndEdg = dicMesh.boundaryEdges ;
        % Element 2 edge connectivity (interior edges only)
            elem2edge = tri2edg ;
            elem2edge(bndEdg,:) = 0 ;
        % Retrieve indices
            [edg,ele] = find(elem2edge) ;
            [edg,is] = sort(edg) ;
            ele = ele(is) ;
        % Change values to [-1 1] for differenciation
            val = repmat([-1;1],[numel(edg)/2 1]) ;
        % Differenciation
            De = sparse(edg,ele,val,dicMesh.nEdges,dicMesh.nElems) ;
        % Normalize with centroid distance (finite difference)
            dC = De*dicMesh.centroid ;
            dC = sqrt(sum(dC.^2,2)) ;
            idC = 1./dC(:) ;
            De(~bndEdg,:) = idC(~bndEdg).*De(~bndEdg,:) ;
        % Integration: use element size
            areas = dicMesh.elemSize ;
            meanArea = sparse(edg,ele,1/2,dicMesh.nEdges,dicMesh.nElems)*areas ;
            De = sqrt(meanArea).*De ;
    % Strain Criterion Matrix
        Ed = blkdiag(De,De,De,De) ; % [4*nEdges 4*nElems]
    
% STRAIN GRADIENT CRITERION: BULK (2nd-order elements)

    % Second gradient
        G2 = dicMesh.grad2Mat ;
    % Integration weights
        for cc = 1:numel(G2) ; G2{cc} = sqrt(areas).*G2{cc} ; end
    % Matrix
        O = sparse(dicMesh.nElems,dicMesh.nNodes) ;
        G2 = [...
                G2{1,1} O ; ... % d2u1_dx1dx1
                G2{2,2} O ; ... % d2u1_dx2dx2
                G2{1,2} O ; ... % d2u1_dx1dx2
                O G2{1,1} ; ... % d2u2_dx1dx1
                O G2{2,2} ; ... % d2u2_dx2dx2
                O G2{1,2} ; ... % d2u2_dx1dx2
             ] ;    
    
    
%% PERFORM DIC !

% PARAMETERS
    % Displacement guess
        startWithNavDICPositions = 'none' ; % Use a preceding computation as guess: 'all', 'none' or a vector of frames
        addPreviousCorrection = true ; % When possible, add the previous correction (velocity or difference with navDIC positions) to the initialization
    % Reference Image 
        weightCurrentImage = 0.05 ; 0.025 ; % After convergence, add the current image to the reference image ([0->1])
    % Image gradient estimation and smoothing
        kernelModel =   ... 'finiteDiff' ... first order finite difference
                         'gaussian' ... optimized gaussian
                        ... 'cos2' ... hamming window
                        ;
        sizeImageKernel = 1 ; % Size of the derivation kernel if needed (allows smoothing)
    % Image Warping
        imWarpInterpOrder = 'linear' ;
    % Image difference criterion
        diffCriterion = ... 'Diff' ... Simple difference
                        ... 'ZM_Diff' ... Zero-mean difference
                         'ZM_N_Diff' ... Normalized Zero-mean difference
                         ;
    % Descent Algorithm
        method = ... 'full-GN' ... full Gauss-Newton
                  'mod-GN' ... modified (assume "grad(g(x+u))=grad(f(x))") OK for small perturbations
                 ... 'quasi-GN' ... not implemented
                 ... 'ICGN' ... not implemented
                 ... 'FCGN' ... not implemented
                    ;
    % Geometry validation criteria
        cullOutOfFrame = true ; % Cull out of frame points
        WEIGHT = INSIDE ; MAPPING ; % % % For local averaging and difference image moments computations
        minCorrCoeff = 0 ; % Below this, elements are culled
        maxMeanElemResidue = Inf ; % Above this, elements are culled
        thresholdValidGeometry = 10 ; % Check correlation. coeffs when the (normA/minNorm)<thresholdValidGeometry. (0 disable the check)
    % Regularization
        stepRatio = 1 ; %0.15 ; % Descent step ratio, damping the convergence
        regCrit = 'rel' ; % second gradient minimization: absolute variation ('abs') or relative ('rel')
        beta = 1*1e7 ; % Strain gradient penalisation coefficient
        epsTrsh = 1e0 ; % Limit value for the regularisation weights (active when regCrit = 'rel')
    % Convergence Criteria
        maxIt = 1000 ; % Maximum number of Newton-Raphson iterations
        minNorm = 1e-1 ; % Maximum displacement of a node
        maxResidueRelativeVariation = Inf ; -.001 ; % Maximum relative variation of the image residue (RMSE)
        minCorrdU = -.999 ; % Maximum update correlation between two consecutive iterations (avoids oscillations)
    % Displacement Processing
        exportTOnavDIC = true ;
        reverseReference = true ;
        strainOnNodes = true ;
    % Plotting
        plotRate = 1 ; % Plot Refresh Frequency 
        plotEachIteration = false ; % Plot at every iteration (without necessary pausing, bypass plotRate)
        plotEachFrame = false ; % Plot at every Frame end (without necessary pausing, bypass plotRate)
        pauseAtPlot = false ; % Pause at each iteration for debugging
    % Watch CPU 
        codeProfile = false ;
    
    
% RUN DIC
    % Initialize
        globalDIC_04_ImageKernels ;
        globalDIC_05_InitDIC ;
    % Run
        if codeProfile ; profile on ; end
        globalDIC_06_PerformDIC ;
        if codeProfile; profile off ; profile off ; end
    % Send to navDIC
        if exportTOnavDIC
            seed.MovingPoints(:,:,frames) = Xn ;
            seed.DataFields = [] ;
        end
