%% CREATE A GRID MESH SEED
    global hd
    cam = 1 ;
    refFrame = hd.CurrentFrame ; 1:hd.nFrames ;
    
    if numel(refFrame)==1
        refImg = hd.Images{cam}{refFrame} ;
    else
        refImg = zeros(size(hd.Images{cam}{end})) ;
        for ii = 1:numel(refFrame)
            refImg = refImg + double(hd.Images{cam}{ii}) ;
        end
        refImg = refImg/numel(refFrame) ;
        refImg = cast(refImg,class(hd.Images{cam}{end})) ;
    end
    [nI,nJ,nC] = size(refImg) ;
        
    clf ;
    imagesc(refImg)
    axis ij
    axis equal
    axis tight
    colormap gray
    
%% REGULAR GRID
    % Macro parameters
        N = [9 9] ; % number of unit cells
        pos0 = round( [1 1 ; nJ 1 ; nJ nI ; 1 nI] + 0.02*[1 1 ; -1 1 ; -1 -1 ; 1 -1].*[nJ nI] ) ; % initial grid corner position
        margin = inf*[1 1 1 1] ; % [left right top bottom]
        elmt =  ...pkg.geometry.mesh.elements.LagrangeElement('quad',2) ...
                pkg.geometry.mesh.elements.base.Quadrangle ...
                ;
    % Micro parameters
        [ee,nn] = meshgrid((0:N(1))/N(1),(0:N(2))/N(2)) ; % unit cell corners
        switch 2
            case 0 % regular grid
                idx = reshape(1:numel(ee),size(ee)) ; % unit cell corner indices
                p1 = idx(1:end-1,1:end-1) ; p2 = p1+1 ; p3 = p2+N(2)+1 ; p4 = p3-1 ;
                elems = [p1(:) p2(:) p3(:) p4(:)] ;
            case 1 % honeycomb lattices
                angle = 60*pi/180 ; margin = [0/10 2/10 0/10 5/10 3/10]+eps ; % [left right top bottom X1+X2]
                angle = 112*pi/180 ; margin = [0/10 0/10 2/10 2/10 3/10]+eps ; % [X1min X1max X2min X2max X1+X2]
                ca = cos(angle)/(1+cos(angle)) ;
                de = .5*[0 ca 1 1+ca] ;
                ee = ee + reshape(de/N(1),1,1,[]) ;
                nn = nn + reshape([0 .5 .5 0]/N(2),1,1,[]) ;
            % Bars
                idx = reshape(1:numel(ee),size(ee)) ; % unit cell corner indices
                elems = [idx(:,:,[1 2]) ; idx(:,:,[2,3]) ; idx(:,:,[3,4])] ;
                elems = reshape(elems,[],2) ;
                elems = [elems ; reshape([idx(:,1:end-1,4) idx(:,2:end,1)],[],2)] ;
                elems = [elems ; reshape([idx(1:end-1,:,2) idx(2:end,:,1)],[],2)] ;
                elems = [elems ; reshape([idx(1:end-1,:,3) idx(2:end,:,4)],[],2)] ;
            case 2 % chiral lattice
                re = 20/100 ;
                ee = ee + cat(3,0,.5+re,.5,.5-re,.5)./N(1) ;
                nn = nn + cat(3,0,.5,.5-re,.5,.5+re)./N(2) ;
                idx = reshape(1:numel(ee),size(ee)) ; % unit cell corner indices
                elems = [...
                            idx(1:end-1,1:end-1,[2 3]) ...
                            idx(1:end-1,1:end-1,[3 4]) ...
                            idx(1:end-1,1:end-1,[4 5]) ...
                            idx(1:end-1,1:end-1,[5 2]) ...
                            idx(1:end-1,1:end-1,[2 4]) ...
                            idx(1:end-1,1:end-1,[3 5]) ...
                            idx(1:end-1,1:end-1,[1 3]) ...
                            cat(3,idx(1:end-1,1:end-1,2),idx(1:end-1,2:end,1)) ...
                            cat(3,idx(1:end-1,1:end-1,4),idx(2:end,1:end-1,1)) ...
                            cat(3,idx(1:end-1,1:end-1,5),idx(2:end,2:end,1)) ...
                        ] ;
        end
    % Cull out of frame
        inX1 = ee(:)>=-margin(1)./N(1) & ee(:)<=1+margin(2)./N(1) ; 
        inX2 = nn(:)>=-margin(3)./N(2) & nn(:)<=1+margin(4)./N(2) ;
        inX1X2 = true ; %(ee(:)*nX+nn(:)*nY-nX-nY)<=margin(5) ;
        valid = inX1 & inX2 & inX1X2 ; 
        elems(any(~valid(elems),2),:) = [] ; % keep bars with all valid points
        valid = valid & ismember(1:numel(ee(:)),elems)' ; % keep points linked to at least one element
        ee = ee(valid) ; nn = nn(valid) ;
        cv = cumsum(valid) ;
        elems = cv(elems) ;
    % Initial grid corners position
        poly = findobj(gca,'type','images.roi.polygon') ;
        if isempty(poly) 
            pos0 = pkg.geometry.mesh.elements.base.Quadrangle().evalAt(elmt.NodeLocalCoordinates)*pos0(1:4,:) ; 
            poly = images.roi.Polygon('Parent',gca ...
                                    ,'Position',pos0 ...
                                    ,'FaceSelectable',false ...
                                    ,'FaceAlpha',0 ...
                                    ) ; 
            poly.Position = pos0 ;
        end
    % Grid warping function
        Nx = elmt.evalAt([ee(:) nn(:)]) ;
        grid = @(pos) Nx*pos ;
    % Point plot
        elems = reshape(elems,[],size(elems,ndims(elems)));
        pts = findobj(gca,'type','patch') ;
        if isempty(pts) 
            pts = patch(gca,'faces',elems ...
                            ,'vertices',grid(poly.Position) ...
                            ,'hittest','off' ...
                            ,'pickableparts','none' ...
                            ,'facecolor','none' ...
                            ,'edgecolor','r' ...
                            ,'marker','.' ...
                            ,'markeredgecolor','r' ...
                            ,'markersize',25 ...
                            ) ; 
        end
        pts.Faces = elems ;
        pts.Vertices = grid(poly.Position) ;
    % Polygon listener
        list = event.listener(poly,'MovingROI',@(src,evt)set(pts,'vertices',grid(poly.Position))) ;
        
%% Adjust points individually
    delete(findobj(gca,'type','images.roi.point'))
    % Create the control points
        nodes = images.roi.Point.empty ;
        for ii = 1:nX*nY ; nodes(end+1) = images.roi.Point(gca,'Position',pts.Vertices(ii,:)) ; end
    % Listen their position
        for ii = 1:nX*nY ; addlistener(nodes(ii),'MovingROI',@(src,evt)set(pts,'vertices',reshape([nodes.Position],[],nX*nY)')) ; end
        
%% Create the seed
    newSeed = copy(hd.Seeds(end)) ;
    newSeed.Name = 'Chiral' ;
    newSeed.Elems = elems ;
    newSeed.Points = pts.Vertices ;
    newSeed.MovingPoints = repmat(pts.Vertices,[1 1 hd.nFrames])*NaN ;
    newSeed.CamIDs = cam ;
    newSeed.refImgs = {refImg} ;
    
%% Push in navDIC
    hd.Seeds(end+1) = newSeed ;
    
%% Backup info
    hd.UserData = [] ;
    hd.UserData.(newSeed.Name).nX = nX ; 
    hd.UserData.(newSeed.Name).nY = nY ; 
    hd.UserData.(newSeed.Name).ee = ee ; 
    hd.UserData.(newSeed.Name).nn = nn ; 
    hd.UserData.(newSeed.Name).pos0 = pos0 ; 
    hd.UserData.(newSeed.Name).polypos = poly.Position ; 
    hd.UserData.(newSeed.Name).angle = angle ; 
    hd.UserData.(newSeed.Name).margin = margin ; 
    hd.UserData.(newSeed.Name).valid = valid ; 
    
    
%% Assemble multiple grids
    seeds = hd.Seeds(3:end) ;
    
    newSeed = copy(seeds(end)) ;
    newSeed.Name = 'AllPoints' ;
    newSeed.Elems = cat(1,seeds.Elems) ;
    newSeed.Points = cat(1,seeds.Points) ;
    newSeed.MovingPoints = cat(1,seeds.MovingPoints) ;
    
    % Delaunay triangulation
    newSeed.Elems = delaunay(newSeed.Points) ;
    xe = reshape(newSeed.Points(newSeed.Elems,:),[],3,2) ;
    newSeed.Elems(polyarea(xe(:,:,1),xe(:,:,2),2)<1,:) = [] ;
    
