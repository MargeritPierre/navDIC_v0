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
    % Number of nodes
        nX = 9 ; round(150/25) ; 
        nY = 9 ; NaN ; round(nI/nJ*(nX-1))+1 ;
    % Initial grid corners position
        pos0 = round( [1 1 ; nJ 1 ; nJ nI ; 1 nI] + 0.02*[1 1 ; -1 1 ; -1 -1 ; 1 -1].*[nJ nI] ) ;
        poly = findobj(gca,'type','images.roi.polygon') ;
        if isempty(poly) 
            poly = images.roi.Polygon('Parent',gca ...
                                    ,'Position',pos0 ...
                                    ,'FaceSelectable',false ...
                                    ,'FaceAlpha',0 ...
                                    ) ; 
        end
    % Grid Function
        if isnan(nY)
            L = range(poly.Position,1) ;
            nY = round(L(2)/L(1)*(nX-1))+1 ;
        end
        [ee,nn] = meshgrid(linspace(0,1,nX),linspace(0,1,nY)) ;
        N = [(1-ee(:)).*(1-nn(:)) ee(:).*(1-nn(:)) ee(:).*nn(:) (1-ee(:)).*nn(:)] ; 
        grid = @(pos) N*pos ;
    % Point plot
        IND = reshape(1:nX*nY,nY,nX) ;
        p1 = IND(1:end-1,1:end-1) ; p2 = p1+1 ; p3 = p2+nY ; p4 = p3-1 ;
        elems = [p1(:) p2(:) p3(:) p4(:)] ;
        pts = findobj(gca,'type','patch') ;
        if isempty(pts) 
            pts = patch(gca,'faces',elems ...
                            ,'vertices',grid(poly.Position) ...
                            ,'hittest','off' ...
                            ,'pickableparts','none' ...
                            ,'facecolor','none' ...
                            ,'edgecolor','b' ...
                            ,'marker','.' ...
                            ,'markeredgecolor','b' ...
                            ,'markersize',25 ...
                            ) ; 
        end
        pts.Faces = elems ;
        pts.Vertices = grid(poly.Position) ;
    % Polygon listener
        list = event.listener(poly,'MovingROI',@(src,evt)set(pts,'vertices',grid(poly.Position))) ;
        
%% HEXAGONAL GRID
    % Number of cells
        nX = 28 ; nY = 18 ; angle = 111*pi/180 ;
        pos0 = round( [1 1 ; nJ 1 ; nJ nI ; 1 nI] + 0.02*[1 1 ; -1 1 ; -1 -1 ; 1 -1].*[nJ nI] ) ;
        elmt = pkg.geometry.mesh.elements.base.Quadrangle ;
    % Presets
        nX = 18 ; nY = 30 ; angle = 60*pi/180 ; pos0 = [422.98 15.412;3114.3 24.466;3104.2 2613.3;415.49 2602.5] ; margin = [0/10 2/10 0/10 5/10 3/10]+eps ; % [X1 X2 X1+X2]
        %nX = 29 ; nY = 19 ; angle = 112*pi/180 ; pos0 = [403.1 82.83;2997 94.98;2984 2659;388.8 2648] ; margin = [0/10 0/10 2/10 2/10 3/10]+eps ; % [X1min X1max X2min X2max X1+X2]
        elmt =  ...pkg.geometry.mesh.elements.LagrangeElement('quad',2) ...
                pkg.geometry.mesh.elements.base.Quadrangle ...
                ;
    % Initial grid corners position
        poly = findobj(gca,'type','images.roi.polygon') ;
        pos0 = pkg.geometry.mesh.elements.base.Quadrangle().evalAt(elmt.NodeLocalCoordinates)*pos0(1:4,:) ; 
        if isempty(poly) 
            poly = images.roi.Polygon('Parent',gca ...
                                    ,'Position',pos0 ...
                                    ,'FaceSelectable',false ...
                                    ,'FaceAlpha',0 ...
                                    ) ; 
        end
        poly.Position = pos0 ;
    % Grid Function
        [ee,nn] = meshgrid((0:nX)/nX,(0:nY)/nY) ; % initial corners
        ca = cos(angle)/(1+cos(angle)) ;
        de = .5*[0 ca 1 1+ca] ;
        ee = ee + reshape(de/nX,1,1,[]) ;
        nn = nn + reshape([0 .5 .5 0]/nY,1,1,[]) ;
    % Bars
        idx = reshape(1:numel(ee),size(ee)) ;
        bars = [idx(:,:,[1 2]) ; idx(:,:,[2,3]) ; idx(:,:,[3,4])] ;
        bars = reshape(bars,[],2) ;
        bars = [bars ; reshape([idx(:,1:end-1,4) idx(:,2:end,1)],[],2)] ;
        bars = [bars ; reshape([idx(1:end-1,:,2) idx(2:end,:,1)],[],2)] ;
        bars = [bars ; reshape([idx(1:end-1,:,3) idx(2:end,:,4)],[],2)] ;
    % Cull out of frame
        inX1 = ee(:)>=-margin(1)./nX & ee(:)<=1+margin(2)./nX ; 
        inX2 = nn(:)>=-margin(3)./nY & nn(:)<=1+margin(4)./nY ;
        inX1X2 = (ee(:)*nX+nn(:)*nY-nX-nY)<=margin(5) ;
        valid = inX1 & inX2 & inX1X2 ; 
        bars(any(~valid(bars),2),:) = [] ; % keep bars with all valid points
        valid = valid & ismember(1:numel(ee(:)),bars)' ; % keep points linked to at least a bar
        ee = ee(valid) ; nn = nn(valid) ;
        cv = cumsum(valid) ;
        bars = cv(bars) ;
    % Grid warping function
        %N = [(1-ee(:)).*(1-nn(:)) ee(:).*(1-nn(:)) ee(:).*nn(:) (1-ee(:)).*nn(:)] ; 
        N = elmt.evalAt([ee(:) nn(:)]) ;
        grid = @(pos) N*pos ;
    % Point plot
        elems = reshape(bars,[],2);%(1:numel(nn))' ;%[p1(:) p2(:) p3(:) p4(:)] ;
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
    newSeed.Name = 'HEX_warp' ;
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
    
