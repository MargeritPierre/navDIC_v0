%% CREATE A NEW SEED CORRESPONDING TO UNIT CELLS
    
%% CREATE NEW

% Parameters
    % New Seed
        nodesSeedNumber = 2 ;
        newSeedNumber = 5 ;
        newSeedName = 'UnitCells2' ;
    % Unit Cell Nodes
        uNodJ = [1:2:9] + 0.5 ;
        uNodI = [2:2:8] + 0.5 ;
    
% Create the new Seed
    nodesSeed = hd.Seeds(nodesSeedNumber) ;
    
% Find the Nodes grid shape (n1,n2)
    Nodes = nodesSeed.MovingPoints(:,:,1) ;
    % First Direction (Approximated with the first difference)
        dx1 = diff(Nodes(1:2,:),1,1) ;
        v1 = dx1/norm(dx1) ;
        L1 = sum((Nodes(end,:)-Nodes(1,:)).*v1) ;
        Ni = round(L1/norm(dx1))+1 ;
    % Second Direction
        Nj = size(Nodes,1)/Ni ;
    
% Nodes Coordinates as Grid
    [jj,ii] = meshgrid(1:Nj,1:Ni) ;
    xx = reshape(Nodes(:,1),[Ni Nj]) ;
    yy = reshape(Nodes(:,2),[Ni Nj]) ;
    
% New nodes coordinates
    uII = uNodI(:) + zeros(size(uNodJ(:)')) ;
    uJJ = zeros(size(uNodI(:))) + uNodJ(:)' ;
    uXX = interp2(jj,ii,xx,uJJ(:),uII(:)) ;
    uYY = interp2(jj,ii,yy,uJJ(:),uII(:)) ;
    uNodes = [uXX(:) uYY(:)] ;
    
% Set the new seed
    newSeed = copy(nodesSeed) ;
    newSeed.Name = newSeedName ;
    hd.Seeds(newSeedNumber) = newSeed ;
    newSeed.Points = uNodes ;
    % Create quad mesh
        nJ = numel(uNodJ) ; nI = numel(uNodI) ;
        newSeed.Elems = (1:nI-1)' + [0 nI nI+1 1] ;
        newSeed.Elems = repmat(newSeed.Elems,[nJ-1 1]) + nI*repelem((0:nJ-2)',nI-1) ;
    
    
% Interpolate to the new nodes
    interpMat = nodesSeed.interpMat(uNodes) ;
    newSeed.MovingPoints = reshape(interpMat*nodesSeed.MovingPoints(:,:),[],2,size(nodesSeed.MovingPoints,3)) ;
    newSeed.computeDataFields ;










        
        
        