function seed2D = Jauge2Dfrom3D( seed3D, hd )

for s = 1:length( seed3D )
    for cam = 1:2 
        ind = (s-1)*2 + cam ;
        seed2D(ind) = navDICSeed_2D_Jauge( hd ) ;
        seed2D( ind).Name = [seed3D(s).Name,'_cam',num2str(cam)] ;
        seed2D(ind).CamIDs = seed3D(s).CamIDs( cam ) ;
        seed2D(ind).refImgs = seed3D(s).refImgs( cam ) ;
        seed2D(ind).RefFrame = seed3D(s).RefFrame ;
        seed2D(ind).Points = seed3D(s).Points( :, :, cam ) ;
        seed2D(ind).MovingPoints = seed3D(s).MovingPoints( :, :, :, cam ) ;
        seed2D(ind).Displacements = seed3D(s).Displacements( :, :, :, cam ) ;
        seed2D(ind).Strains = seed3D(s).Strains( :, 1, :, cam ) ;
        seed2D(ind).corrSize = seed3D(s).corrSize( cam ) ;
        seed2D(ind).L0 = seed3D(s).L0( :, cam ) ;
    end
end