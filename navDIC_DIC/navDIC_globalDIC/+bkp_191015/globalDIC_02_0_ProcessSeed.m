%% GLOBAL DIC: PROCESS THE SEED

% Get the Seed
    Seed = hd.Seeds(seedNumber) ;

% Get Infos
    Elems = Seed.Triangles ;
    nNodesByElements = size(Elems,2) ;
    switch refConfig
        case 'Nodes'
            Nodes = Seed.Points ;
        case 'Current'
            Nodes = Seed.MovingPoints(:,:,frames(refFrame)) ;
            if any(isnan(Nodes(:))) ; error('The Current Node Configuration is not defined') ; end
    end
    nElems = size(Elems,1) ;
    nNodes = size(Nodes,1) ;
    
% Processes
    globalDIC_02_1_ProcessGeometries ;
    globalDIC_02_2_ShapeFunctions ;
    
    