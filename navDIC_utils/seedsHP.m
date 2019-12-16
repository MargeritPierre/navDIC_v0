function seedHP = seedsHP(seed, hd, ddo, om ) 

if isempty( seed ); return; end
if isempty(seed.Uhp) ; return; end 

nbImgs = size(seed.MovingPoints,3) ;
seedHP.CamIDs = seed.CamIDs ;
seedHP.Name = seed.Name ;
seedHP.MovingPoints = zeros( size(seed.MovingPoints) ) ;
Cam = hd.Cameras(seed.CamIDs).Properties ;
try 
    Cam.do = Cam.do(1) + ddo ;
    Cam.pixObjratio = Cam.fpix(1) / Cam.do ;
end

for i = 1:nbImgs
    theta = seed.theta(i) + om ;
    seedHP.MovingPoints(:,:,i) = homoRotTr( seed.MovingPoints(:,:,i), seed.IDrefHP,...
        Cam, theta, seed.Uhp(i) ) ; 
end

seedHP.defUhp = seed.Uhp / Cam.do ;
seedHP.theta = seed.theta ;
seedHP.defTheta = cos( seed.theta ) / cos( seed.theta(1) ) - 1 ;

seedHP.Points = seedHP.MovingPoints(:,:,1) ;
seedHP.Displacements = seedHP.MovingPoints - repmat( seedHP.MovingPoints(:,:,1), [ 1, 1, nbImgs ] ) ;
seedHP.L0 = sum( ( seedHP.MovingPoints( 2:2:end, :, 1) - seedHP.MovingPoints( 1:2:end, :, 1 ) ).^2, 2 ).^.5 ;
L = sum( ( seedHP.MovingPoints(2:2:end,:,: ) - seedHP.MovingPoints(1:2:end,:, :) ).^2, 2 ).^.5 ;
seedHP.Strains  = ( L - repmat( seedHP.L0, [ 1, 1, nbImgs ] ) ) ./ repmat( seedHP.L0, [ 1, 1, nbImgs ] ) ;
