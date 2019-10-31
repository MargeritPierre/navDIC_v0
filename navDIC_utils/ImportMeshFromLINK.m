%% IMPORT A MESH CREATED WITH LINK IN AN EXISTING NAVDIC SESSION
% CREATE A NEW SEED WITH THE MESH
% /!\ FIRST CENTER-CLICK THE LINK NODE CONTAINING THE MESH !
% THEN type "mesh = ans.Data"

newSeedName = 'From_LINK_EdgLen20' ;
seedModelNumber = 5 ;

global hd
% Initialize with a navDIC seed
    newSeed = hd.Seeds(seedModelNumber) ;
    newSeed.Name = newSeedName ;
% Fill the new Seed
    newSeed.Triangles = mesh.Faces ;
    newSeed.Points = mesh.Vertices ;
% Initialize the seed computations
    nP = size(mesh.Vertices,1) ;
    newSeed.MovingPoints = NaN*ones(nP,2,hd.nFrames) ;
    newSeed.Displacements = NaN*ones(nP,2,hd.nFrames) ;
    newSeed.Strains = NaN*ones(nP,3,hd.nFrames) ;
    newSeed.MajorStrains = NaN*ones(nP,1,hd.nFrames) ;
    newSeed.MinorStrains = NaN*ones(nP,1,hd.nFrames) ;
    newSeed.MaxShear = NaN*ones(nP,1,hd.nFrames) ;
    newSeed.PrincipalAngle = NaN*ones(nP,1,hd.nFrames) ;

% Display
    newSeed
    
%% PUSH TO navDIC setup
    hd.Seeds(end+1) = newSeed ;