%% PROJECT DISPLACEMENTS COMPUTED ON A MESH To ANOTHER MESH


SourceSeedNumber = 5 ;
DestinationSeedNumber = 4 ;

global hd

source = hd.Seeds(SourceSeedNumber) ;
destination = hd.Seeds(DestinationSeedNumber) ;

T = source.interpMat(destination.Points) ;


hd.Seeds(DestinationSeedNumber).MovingPoints = reshape(T*reshape(source.MovingPoints,size(source.Points,1),[]),size(destination.Points,1),2,hd.nFrames) ;
hd.Seeds(DestinationSeedNumber).Displacements = reshape(T*reshape(source.Displacements,size(source.Points,1),[]),size(destination.Points,1),2,hd.nFrames) ;
hd.Seeds(DestinationSeedNumber).Strains = ones(size(destination.Points,1),3,hd.nFrames)*NaN ; % Strains will need to be re-computed...
hd.Seeds(DestinationSeedNumber).computeDataFields() ;