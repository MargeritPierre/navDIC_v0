%% PROJECT DISPLACEMENTS COMPUTED ON A MESH TO ANOTHER MESH


SourceSeedNumber = 5 ;
DestinationSeedNumber = 4 ;

global hd

srcSeed = hd.Seeds(SourceSeedNumber) ;
dstSeed = hd.Seeds(DestinationSeedNumber) ;

T = srcSeed.interpMat(dstSeed.Points) ;


hd.Seeds(DestinationSeedNumber).MovingPoints = reshape(T*reshape(srcSeed.MovingPoints,size(srcSeed.Points,1),[]),size(dstSeed.Points,1),2,hd.nFrames) ;
hd.Seeds(DestinationSeedNumber).Displacements = reshape(T*reshape(srcSeed.Displacements,size(srcSeed.Points,1),[]),size(dstSeed.Points,1),2,hd.nFrames) ;
hd.Seeds(DestinationSeedNumber).Strains = ones(size(dstSeed.Points,1),3,hd.nFrames)*NaN ; % Strains will need to be re-computed...
hd.Seeds(DestinationSeedNumber).computeDataFields() ;



%% NEW VERSION (PKG.GEOMETRY.MESH REQUIRED !)

SourceSeedNumber = 3 ;
DestinationSeedNumber = 4 ;
extrap = true ;

global hd

srcSeed = hd.Seeds(SourceSeedNumber) ;
dstSeed = hd.Seeds(DestinationSeedNumber) ;

srcMesh = pkg.geometry.mesh.Mesh('Nodes',srcSeed.Points,'Elems',srcSeed.Elems) ;

%interpMat(this,P,ie,features,extrap,X,tol)
M = srcMesh.interpMat(dstSeed.Points,[],srcMesh.Elems,extrap) ;
[nPtsSrc,nCoord,nFr] = size(srcSeed.MovingPoints) ;
dstSeed.MovingPoints = reshape(M*srcSeed.MovingPoints(:,:),[],nCoord,nFr) ;
dstSeed.computeDataFields ;



