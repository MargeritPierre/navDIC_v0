if 1 % USE THIS TO GO DIRECTLY TO DIC

    % CLEAN WORKSPACE
        clc
        global hd
        clearvars -except hd

    % INITIALIZATION PARAMETERS
        camID = 1 ;
        seedNumber = 1 ;
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
        addPreviousVelocity = true ; % When possible and no navDIC results available (or not used), add the previous motion as convergence help
    % Reference Image 
        weightCurrentImage = 0.3 ; % After convergence, add the current image to the reference image ([0->1])
    % Image gradient estimation and smoothing
        kernelModel =   ... 'finiteDiff' ... first order finite difference
                         'gaussian' ... optimized gaussian
                        ... 'cos2' ... hamming window
                        ;
        sizeImageKernel = 3 ; % Size of the derivation kernel if needed (allows smoothing)
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
        maxMeanElemResidue = Inf ; % Above this, elements are culled
        thresholdValidGeometry = 0 ; % Check correlation. coeffs when the (normA/minNorm)<thresholdValidGeometry. (0 disable the check)
    % Regularization
        beta = 1*1e2 ; % Strain gradient penalisation coefficient
    % Convergence Criteria
        maxIt = 100 ; % Maximum number of Newton-Raphson iterations
        minNorm = 1e-3 ; % Maximum displacement of a node
        minCorrdU = -.9999 ; % Maximum update correlation between two consecutive iterations (avoids oscillations)
    % Displacement Processing
        exportTOnavDIC = false ;
        reverseReference = true ;
        strainOnNodes = false ;
    % Plotting
        plotRate = 1 ; % Plot Refresh Frequency 
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
    % Process
        globalDIC_07_AfterDIC ;
    % Send to navDIC
        if exportTOnavDIC ; globalDIC_08_TOnavDIC ; end

    