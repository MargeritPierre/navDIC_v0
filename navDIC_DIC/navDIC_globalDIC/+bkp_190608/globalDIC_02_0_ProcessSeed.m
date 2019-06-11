%% GLOBAL DIC: PROCESS THE SEED

% Get the Seed
    Seed = hd.Seeds(seedNumber) ;

% Get Infos
    Elems = Seed.Triangles ;
    nNodesByElements = size(Elems,2) ;
    Nodes = Seed.Points ;
    nElems = size(Elems,1) ;
    nNodes = size(Nodes,1) ;
    
% Processes
    globalDIC_02_1_ProcessGeometries ;
    globalDIC_02_2_ShapeFunctions ;
    
    