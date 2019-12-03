if 1 % USE THIS TO GO DIRECTLY TO DIC

    % CLEAN WORKSPACE
        clc
        global hd
        clearvars -except hd

    % INITIALIZATION PARAMETERS
        camID = 1 ;
        seedNumber = 3 ;
        frames = '[1:end]' ; % Frames taken for DIC (allows decimation)
        dicDir = 1 ; % DIC running direction ('forward=1' or 'backward=-1')
        refFrame = 'first' ; % Reference image ('first' , 'last' or number)
        refConfig = 'Nodes' ; % Reference configuration: 'Nodes' (as meshed) or 'Current' (uses preceding computed displacement)
        averagePreviousFrames = true ; % Ref frame is the average of the previous/next ones in forward/backward modes
        normToImageClassRange = true ; % Normalize images to their dataclass range
        timeMeanLength = 0 ; % Time averaging of images
        strainCriterion = 'full' ; % strain gradient penalization: 'full' or 'normal'
        showInit = false ;
        codeProfile = false ; % Code timing
        figTag = 'Global DIC' ;

    % Perform Initialization
        globalDIC_01_LoadFrames ;
        globalDIC_02_0_ProcessSeed ;
        globalDIC_03_0_StrainCriterion
        if showInit ; globalDIC_03_1_ShowInitialization ; end
        
end % END OF INITIALIZATION
    
%% PERFORM DIC !

% PARAMETERS
    % Displacement guess
        startWithNavDICPositions = 'all' ; % Use a preceding computation as guess: 'all', 'none' or a vector of frames
        addPreviousCorrection = true ; % When possible, add the previous correction (velocity or difference with navDIC positions) to the initialization
    % Reference Image 
        weightCurrentImage = 0.025 ; %0.025 ; % After convergence, add the current image to the reference image ([0->1])
    % Image gradient estimation and smoothing
        kernelModel =   ... 'finiteDiff' ... first order finite difference
                         'gaussian' ... optimized gaussian
                        ... 'cos2' ... hamming window
                        ;
        sizeImageKernel = 3 ; % Size of the derivation kernel if needed (allows smoothing)
    % Image Warping
        imWarpInterpOrder = 'cubic' ;
    % Image difference criterion
        diffCriterion =  'Diff' ... Simple difference
                        ... 'ZM_Diff' ... Zero-mean difference
                        ... 'ZM_N_Diff' ... Normalized Zero-mean difference
                         ;
    % Descent Algorithm
        method =  'full-GN' ... full Gauss-Newton
                 ... 'mod-GN' ... modified (assume "grad(g(x+u))=grad(f(x))") OK for small perturbations
                 ... 'quasi-GN' ... not implemented
                 ... 'ICGN' ... not implemented
                 ... 'FCGN' ... not implemented
                    ;
    % Geometry validation criteria
        cullOutOfFrame = true ; % Cull out of frame points
        WEIGHT = INSIDE ; % MAPPING ; % For local averaging and difference image moments computations
        minCorrCoeff = .0 ; % Below this, elements are culled
        maxMeanElemResidue = Inf ; % Above this, elements are culled
        thresholdValidGeometry = 0 ; % Check correlation. coeffs when the (normA/minNorm)<thresholdValidGeometry. (0 disable the check)
    % Regularization
        stepRatio = 1 ;%0.15 ; % Descent step ratio, damping the convergence
        regCrit = 'rel' ; % second gradient minimization: absolute variation ('abs') or relative ('rel')
        beta = 1*1e3 ; % Strain gradient penalisation coefficient
        epsTrsh = 1e0 ; % Limit value for the regularisation weights (active when regCrit = 'rel')
    % Convergence Criteria
        maxIt = 100 ; % Maximum number of Newton-Raphson iterations
        minNorm = 1e-3 ; % Maximum displacement of a node
        maxResidueRelativeVariation = -.001 ; % Maximum relative vriation of the image residue (RMSE)
        minCorrdU = -.999 ; % Maximum update correlation between two consecutive iterations (avoids oscillations)
    % Displacement Processing
        exportTOnavDIC = true ;
        reverseReference = true ;
        strainOnNodes = true ;
    % Plotting
        plotRate = 0 ; % Plot Refresh Frequency 
        plotEachIteration = false ; % Plot at every iteration (without necessary pausing, bypass plotRate)
        plotEachFrame = true ; % Plot at every Frame end (without necessary pausing, bypass plotRate)
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
        if codeProfile; profile viewer ; profile off ; end
    % Send to navDIC
        if exportTOnavDIC
            globalDIC_07_AfterDIC ; 
            globalDIC_08_TOnavDIC ; 
        end

    