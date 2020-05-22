function [in,on] = inpolygons(pts,polys)

    pts = pts(:,:) ; % [nPts nCoord]
    polys = polys(:,:,:) ; % [nPolys nMaxNodes nCoord]
    
    nPts = size(pts,1) ;
    nCoord = size(pts,2) ;
    nPolys = size(polys,1) ;
    nMaxNodes = size(polys,3) ;
    
    polys = permute(polys,[4 3 1 2]) ; % [1 nCoord nPolys nMaxNodes]
    
    % Element Bounding Box (restraint the computations)
        Xmin = min(polys,[],4) ;
        Xmax = max(polys,[],4) ;
        in = all(pts>=Xmin & pts<=Xmax,2)  ;
        
    % Polygons with at least one point inside the BBox
        inPi = any(in,1) ;
        %inPolys = polys(:,:,inPi,:) ; 
        inPolys = polys(:,:,:,:) ; 
        
    % Edges
        edges = cat(5,inPolys,circshift(inPolys,-1,4)) ; % [1 nCoord nPolys nMaxNodes 2]
        
    % Parameter t corresponding to the cross between an edge and an
    % horizontal line starting from each point
        tc = (pts-edges(:,:,:,:,1)).*(1./diff(edges,1,5)) ; % [nPts nCoord nPolys nMaxNodes]
        pc = edges(:,:,:,:,1).*(1-tc) + edges(:,:,:,:,2).*tc ;
        
    % Check % [nPts 1 nPolys nMaxNodes]
        crossY0Right = tc(:,2,:,:,:)>=0 & tc(:,2,:,:,:)<=1 & pc(:,1,:,:,:) >= pts(:,1) ;
        crossY0Left = tc(:,2,:,:,:)>=0 & tc(:,2,:,:,:)<=1 & pc(:,1,:,:,:) <= pts(:,1) ;
        crossX0Top = tc(:,1,:,:,:)>=0 & tc(:,1,:,:,:)<=1 & pc(:,2,:,:,:) >= pts(:,2) ;
        crossX0Bottom = tc(:,1,:,:,:)>=0 & tc(:,1,:,:,:)<=1 & pc(:,2,:,:,:) <= pts(:,2) ;
        
        
    in = reshape(in,[nPts nPolys]) ;    
    on = in ;

return


%% TEST

    pts = [0 0 ; -1 -1 ; -1 1] ;
    polys = [0 0 ; 1 0 ; 1 1 ; 0 1] ;
    polys(:,:,end+1) = [0 0 ; 1 0 ; 1 1 ; 0 1]-2 ;
    polys(:,:,end+1) = [0 0 ; 1 0 ; 1 1 ; NaN NaN] ;
    polys = permute(polys,[3 1 2]) ;

    in = inpolygon(pts,polys)

