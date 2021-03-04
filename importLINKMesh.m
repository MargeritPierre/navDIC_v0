%% IMPORT A MESH CREATED IN LINK

    newSeedNumber = 5 ;
    newSeedName = 'Micro' ;
    meshPort = ans ;

    global hd
    hd.Seeds(newSeedNumber) = copy(hd.Seeds(1)) ;
    newSeed = hd.Seeds(newSeedNumber) ;
    newSeed.Name = newSeedName ;
    
    mesh = meshPort.Data ;
    newSeed.Elems = mesh.Faces ;
    newSeed.Points = mesh.Vertices(:,1:2) ;
    newSeed.MovingPoints = repmat(mesh.Vertices(:,1:2),[1 1 hd.nFrames])*NaN ;
    newSeed.computeDataFields ;
    