classdef navDICSeed_2D_DistMesh < navDICSeed_2D_Surface
   
    properties
        Triangles = [] ;
        Quadrangles = [] ;
        ROI = [] ;
        Shapes = [] ;
        DataOnNodes = false ;
    end
    
    properties
        Elems
    end
    
    properties (Transient,Dependent)
        %Triangles, Quadrangles % <TODO> Later
        Edges
    end
    
    properties (Transient)
        DataFields = struct() ;
    end
    
    
    
methods
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% CONSTRUCT / MODIFY / DELETE THE SEED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function obj = navDICSeed_2D_DistMesh(hd)
        % Initialize
            obj = obj@navDICSeed_2D_Surface(hd) ;
            obj.Class = 'navDICSeed_2D_DistMesh' ;
        % Draw the Shapes defining the region to mesh
            uiContextMenu = {'h0',num2cell(num2str(round((1.4).^([1:17,Inf])')),2),{'29'};...
                'h',num2cell(num2str(round((1.4).^(1:17)')),2),{'29'};...
                'l',num2cell(num2str(round((1.4).^(8:23)')),2),{'79'};...
                } ;
            obj.drawToolH = drawingTool('drawROI',true ...
                                        ,'background',obj.refImgs{1} ...
                                        ,'uicontextmenu',uiContextMenu ...
                                        ,'updateCallback',@(H)navDIC_processShapesForDistMesh(H)...
                                        ) ;
            obj.Points = obj.drawToolH.DistMesh.Points ;
            obj.Triangles = obj.drawToolH.DistMesh.Triangles ;
        % INITIALIZE
            obj.MovingPoints = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
    end


    function obj = modify(obj,hd)
        % ALLOW THE USER TO MODIFY THE MESH
            drawToolH = drawingTool(obj.drawToolH) ;
        % Has the mesh been changed ?
            Points = drawToolH.DistMesh.Points ;
            Tris = drawToolH.DistMesh.Triangles ;
            if size(obj.Points,1)==size(Points,1) && sqrt(sum((obj.Points(:)-Points(:)).^2))<1
                return ; % No need for modifications
            end
        % MODIFY THE MESH DATA
            if all(isnan(obj.MovingPoints(:))) 
                MovingPoints = ones(size(Points,1),1).*obj.MovingPoints(1,:,:) ;
            else
                % ASK THE USER WHAT TO DO
                    answer = questdlg(...
                                {'THE MESH HAS BEEN MODIFIED.',...
                                'WHAT DO YOU WANT TO DO WITH THE EXISTING DATA ?'}...
                                ,'/!\ WARNING'...
                                ,'Project','Clear','Cancel','Project'...
                                ) ; 
                    if strcmp(answer,'Cancel') ; return ; end
                    % Init with NaNs
                        MovingPoints = ones(size(Points,1),2,hd.nFrames)*NaN ;
                    % Project previous computations if wanted
                        if strcmp(answer,'Project')
                            % Interpolation Matrix
                                T = obj.interpMat(Points) ;
                            % Projection of displacements
                                MovingPoints = reshape(T*reshape(obj.MovingPoints,size(obj.Points,1),[]),size(Points,1),2,hd.nFrames) ;
                        end
            end
            % Crush the old mesh
                obj.drawToolH = drawToolH ;
                obj.Points = Points ;
                obj.Triangles = Tris ;
                obj.MovingPoints = MovingPoints ;
                obj.DataFields = [] ;
    end

    function delete(obj)
    % Delete the seed
    end
        
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH-SPECIFIC FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------- MESH GEOMETRY / CONNECTIVITY ---------------------------------

    function set.Elems(obj,elems)
    % Set mesh elements (triangles, quadrangles)
        elems = round(unique(elems,'rows','stable')) ;
        nNotNan = sum(~isnan(elems),2) ;
        obj.Triangles = elems(nNotNan==3,:) ;
        obj.Quadrangles = elems(nNotNan==4,:) ;
        if ~isempty(obj.Triangles) ; obj.Triangles = obj.Triangles(:,1:3) ; 
        else ; obj.Triangles = [] ;
        end
        if ~isempty(obj.Quadrangles) ; obj.Quadrangles = obj.Quadrangles(:,1:4) ;  
        else ; obj.Quadrangles = [] ;
        end
        obj.Elems = elems ;
    end

    function elems = get.Elems(obj)
    % Return all elements of the mesh
        elems = obj.Elems ;
        if ~isempty(obj.Triangles) ; elems = [elems ; obj.Triangles NaN(size(obj.Triangles,1),size(elems,2)-3)] ; end
        if ~isempty(obj.Quadrangles) ; elems = [elems ; obj.Quadrangles NaN(size(obj.Quadrangles,1),size(elems,2)-4)] ; end
        elems = unique(elems,'rows','stable') ;
    end

    function set.Edges(obj,edges)
    % Set the object's edges
        %obj.Edges = round(unique(edges,'rows')) ;
    end

    function edges = get.Edges(obj)
    % Get the object's unique edges
        [edges,~] = getEdges(obj) ;
    end

    function [edges,elem2edg] = getEdges(obj)
    % Get object's edges connectivity
        % Edge nodes
            e1 = obj.Elems ;
            Nn = sum(~isnan(e1),2) ;
            e2 = circshift(e1,-1,2) ;
            e2(sub2ind(size(e2),1:size(e2,1),Nn')) = e2(:,end) ;
        % Sorting 
            edges = [e1(:) e2(:)] ;
            edges = sort(edges,2) ;
        % Corresponding element
            elmt = repmat((1:size(e1,1))',[size(e1,2) 1]) ;
        % Cull invalid edges
            nans = any(isnan(edges),2) ;
            edges = edges(~nans,:) ;
            elmt = elmt(~nans) ;
        % Unique list
            [edges,~,ie] = unique(edges,'rows') ;
            elem2edg = sparse(ie(:),elmt(:),1,size(edges,1),size(e1,1)) ;
    end

    function outE = boundaryEdges(obj)
    % Return a logical array with true for edges on mesh boundaries
        [~,elem2edg] = obj.getEdges ;
        outE = sum(elem2edg,2)==1 ;
    end

    function M = elem2nod(obj,values)
    % Return a matrix so that f = M*F, where f is defined on Nodes and
    % F is defined on elements; size(M) = [nNodes nElems]
        if nargin<2 ; values = 'mean' ; end
        list = obj.Elems ;
        ii = list(:) ;
        jj = repmat([1:size(list,1)]',[size(list,2) 1]) ;
        nans = isnan(ii) | isnan(jj) ;
        M = sparse(ii(~nans),jj(~nans),1,size(obj.Points,1),size(list,1)) ;
        switch values
            case 'ones' % Do nothing
            case 'logical'
                M = logical(M) ;
            case 'mean'
                M = (1./sum(M,2)).*M ;
            otherwise % double values
        end
    end
    

% ------- GLOBAL/LOCAL COORDINATES MAPPING ------------------------------
% Works for triangles and quadrangles only (P1)
% (e1,e2) in [0->1]²
%               3        |          4--3
%        TRI:   | \      |  QUAD:   |  |
%               1--2     |          1--2

    function [nNodes,isTri,isQuad] = getElemTypes(obj)
    % Return the number of valid nodes
        nNodes = sum(~isnan(obj.Elems),2) ;
        isTri = nNodes==3 ;
        isQuad = nNodes==4 ;
    end
    
    function [elmtData,isTri,isQuad] = prepareElements(obj,elmt,nodeData)
    % Prepare the element lists for mapping functions
        if nargin<2 ; elmt = 1:size(obj.Elems,1) ; end
        if nargin<3 ; nodeData = obj.Points ; end
        nodeData = nodeData(:,:) ; nD = size(nodeData,2) ;
    % Data associated to each element node
        elems = obj.Elems(elmt,:) ;
        valid = ~isnan(elems) ; elems(~valid) = 1 ;
        elmtData = reshape(nodeData(elems(:),:),[size(elems) nD]) ;
        elmtData(repmat(~valid(:),[nD 1])) = NaN ;
        elmtData = permute(elmtData,[1 3 2]) ; % [elmt coord node]
    % Element types
        [~,isTri,isQuad] = getElemTypes(obj) ;
        isTri = isTri(elmt) ;
        isQuad = isQuad(elmt) ;
    end

    function X = globalCoordinates(obj,E,elmt,refPts)
    % Return the global coordinates X = [x1 x2] of ...
    % points given by local coordinates E = [e1 e2] ...
    % in elements elmt, with reference nodes refPts
        if nargin<3 ; elmt = 1:size(obj.Elems,1) ; end
        if nargin<4 ; refPts = obj.Points ; end
        elmt = elmt(:) ; E = E(:,:) ; refPts = refPts(:,:) ;
        elmt = elmt + E(:,1)*0 ; E = E + elmt(:,1)*0 ;
    % Prepare elements
        [elmtPts,isTri,isQuad] = prepareElements(obj,elmt,refPts) ;
    % Mapping
        X = NaN(size(E)) ;
        if any(isTri)
            X(isTri,:) =  elmtPts(isTri,:,1).*(1-E(isTri,1)-E(isTri,2)) ...
                            + elmtPts(isTri,:,2).*E(isTri,1) ...
                            + elmtPts(isTri,:,3).*E(isTri,2) ;
        end
        if any(isQuad)
            X(isQuad,:) =  elmtPts(isQuad,:,1).*(1-E(isQuad,1)).*(1-E(isQuad,2)) ...
                            + elmtPts(isQuad,:,2).*E(isQuad,1).*(1-E(isQuad,2)) ...
                            + elmtPts(isQuad,:,3).*E(isQuad,1).*E(isQuad,2) ...
                            + elmtPts(isQuad,:,4).*(1-E(isQuad,1)).*E(isQuad,2) ;
        end
    end

    function [J,detJ,invJ] = jacobian(obj,E,elmt,refPts)
    % Return the jacobian  J = [dx1_de1 dx2_de1 dx1_de2 dx2_de2] ...
    % of the transformation at given local coordinates E = [e1 e2] ...
    % on elements elmt, with reference nodes refPts
        if nargin<3 ; elmt = 1:size(obj.Elems,1) ; end
        if nargin<4 ; refPts = obj.Points ; end
        elmt = elmt(:) ; E = E(:,:) ; refPts = refPts(:,:) ;
        elmt = elmt + E(:,1)*0 ; E = E + elmt(:,1)*0 ;
    % Prepare elements
        [elmtPts,isTri,isQuad] = prepareElements(obj,elmt,refPts) ;
    % Jacobian
        J = NaN(size(elmt,1),4) ;
        if any(isTri) 
            J(isTri,[1 2]) =  elmtPts(isTri,:,2) - elmtPts(isTri,:,2) ; 
            J(isTri,[3 4]) =  elmtPts(isTri,:,3) - elmtPts(isTri,:,2) ; 
        end
        if any(isQuad)
            J(isQuad,[1 2]) =  (1-E(isQuad,2)).*(elmtPts(isQuad,:,2)-elmtPts(isQuad,:,1)) ...
                                + E(isQuad,2).*(elmtPts(isQuad,:,3)-elmtPts(isQuad,:,4)) ;
            J(isQuad,[3 4]) =  (1-E(isQuad,1)).*(elmtPts(isQuad,:,4)-elmtPts(isQuad,:,1)) ...
                                + E(isQuad,1).*(elmtPts(isQuad,:,3)-elmtPts(isQuad,:,2)) ;
        end
    % Determinant
        if nargout>1 ; detJ = J(:,1).*J(:,4)-J(:,2).*J(:,3) ; end
    % Inverse Jacobian
        if nargout>2 ; invJ = (1./detJ).*[J(:,4) -J(:,2) -J(:,3) J(:,1)] ; end
    end

    function H = hessian(obj,E,elmt,refPts)
    % Return the hessian ...
    % H = [d2x1_de1de1 d2x2_de1de1 d2x1_de2de2 d2x2_de2de2 d2x1_de1de2 d2x2_de1de2] ...
    % of the transformation at given local coordinates E = [e1 e2] ...
    % on elements elmt, with reference nodes refPts
        if nargin<3 ; elmt = 1:size(obj.Elems,1) ; end
        if nargin<4 ; refPts = obj.Points ; end
        elmt = elmt(:) ; E = E(:,:) ; refPts = refPts(:,:) ;
        elmt = elmt + E(:,1)*0 ; E = E + elmt(:,1)*0 ;
    % Prepare elements
        [elmtPts,isTri,isQuad] = prepareElements(obj,elmt,refPts) ;
    % Hessian
        H = zeros(size(elmt,1),6) ;
        %if any(isTri) ; H(isTri,:) = 0 ; end
        if any(isQuad)
            %H(isQuad,[1 2 3 4]) =  0 ;
            H(isQuad,[5 6]) =  elmtPts(isQuad,:,1) + elmtPts(isQuad,:,3) - elmtPts(isQuad,:,4) - elmtPts(isQuad,:,2) ;
        end
    end

    function elmt = inside(obj,pts,refPts)
    % Return element indexes in which the given points are included
    % Return NaN if the point is outside all elements
        if nargin<3 ; refPts = obj.Points ; end
        pts = pts(:,:) ;
        elmt = NaN(size(pts,1),1) ; 
    % Edge List
        elems = obj.Elems ;
        edges = cat(3,elems,circshift(elems,1,2)) ;
        ee = repmat((1:size(elems,1))',[1 size(elems,2)]) ;
        validEdg = ~any(isnan(edges),3) ;
    % Global Bounding box
        outside = pts(:,1)>=min(refPts(:,1)) ...
                    & pts(:,1)<=max(refPts(:,1)) ...
                    & pts(:,2)>=min(refPts(:,2)) ...
                    & pts(:,2)<=max(refPts(:,2)) ;
    % Element Bounding Box
        elmt = (1:size(obj.Elems))' ;
        elmtPts = prepareElements(obj,elmt,refPts) ;
        inside = pts(~outside,1)>=min(elmtPts(:,1)) ...
                    & pts(~outside,1)<=max(elmtPts(:,1)) ...
                    & pts(~outside,2)>=min(elmtPts(:,2)) ...
                    & pts(~outside,2)<=max(elmtPts(:,2)) ;
        
        [e1,e2] = obj.localCoordinates(1) ;
        [x1,x2] = obj.globalCoordinates(e1,e2) ;
        x1 = x1([1 3 4 2]) ; x2 = x2([1 3 4 2]) ;
        % In bounding box
            in = pts(:,1)>=min(x1) ...
                    & pts(:,1)<=max(x1) ...
                    & pts(:,2)>=min(x2) ...
                    & pts(:,2)<=max(x2) ;
        % In polygon
            if any(in)
                in(in) = in(in) & inpolygons(pts(in,1),pts(in,2),x1(:),x2(:)) ;
            end
    end

    function [e1,e2] = localize(this,pts,tol,maxIt)
    % Return local coordinates associated to points
        if nargin<3 ; tol = eps*1000*norm(range(this.Position,1)) ; end
        if nargin<4 ; maxIt = 10 ; end
        e1 = NaN(size(pts,1),1) ; e2 = NaN(size(pts,1),1) ;
        % Process only points in the element
            in = this.inside(pts) ;
        % Inverse mapping 
        % see Silva & al, "Exact & efficient interpolation using FE shape functions"
            e1(in) = 0.5 ; e2(in) = 0.5 ;
            for it = 1:maxIt
                % Update candidate coordinates
                    [x1,x2] = this.globalCoordinates(e1(in),e2(in)) ;
                % Residues
                    dx1 = x1-pts(in,1) ;
                    dx2 = x2-pts(in,2) ;
                % Jacobian
                    [dx1_de1,dx2_de1,dx1_de2,dx2_de2] = jacobian(this,e1(in),e2(in)) ;
                % Inverse Jacobian
                    detJ = dx1_de1.*dx2_de2-dx1_de2.*dx2_de1 ;
                    de1_dx1 = dx2_de2./detJ ;
                    de2_dx1 = -dx2_de1./detJ ;
                    de1_dx2 = -dx1_de2./detJ ;
                    de2_dx2 = dx1_de1./detJ ;
                % Increment
                    de1 = - de1_dx1.*dx1 - de1_dx2.*dx2 ;
                    de2 = - de2_dx1.*dx1 - de2_dx2.*dx2 ;
                % Update local coordinates
                    e1(in) = e1(in) + de1 ;
                    e2(in) = e2(in) + de2 ;
                % Increment norm
                    de = sqrt(de1.^2 + de2.^2) ;
                    if de<tol ; break ; end
            end
    end

    function [data,ii,jj] = evalAt(this,pts)
    % Evaluate the pixel values at given points (in global coordinates)
        data = NaN(size(pts,1),size(this.CData,3)) ;
        [e1,e2] = this.localize(pts) ;
        valid = ~isnan(e1) & ~isnan(e2) ;
        [ii,jj] = this.pixelIndices(e1,e2) ;
        ind = sub2ind(size(this.CData),ii(valid).*[1 1 1],jj(valid).*[1 1 1],repmat(1:size(this.CData,3),[size(ii(valid),1) 1])) ;
        data(valid,:) = this.CData(ind) ;
    end
    

% ------- DATA INTERPOLATION -------------------------------------------

    function T = interpMat(obj,Points,extrap,refPoints)
    % Return the interpolation matrix of the mesh on query points so
    % that U(pts) = T*U(Nodes)
        if nargin<3 ; extrap = 'extrap' ; end
        if nargin<4 ; refPoints = obj.Points ; end
        if size(obj.Elems,2)<3 % NO PROPER MESH FOR INTERPOLATION
            % Distance between points and refpoints
                dist = sum((reshape(full(Points),[],1,2)-reshape(full(refPoints),1,[],2)).^2,3) ;
            % Nearest neightbors
                [pt,refPt] = find(dist==min(dist,[],2)) ;
            % Transfer matrix
                T = sparse(pt(:),refPt(:),1,size(Points,1),size(refPoints,1)) ;
        else % THERE IS A MESH FOR INTERPOLATION
            % Convert all elements to triangles (not ideal..)
                elems = obj.Elems ;
                if size(elems,2)>3
                    nDiv = size(elems,2)-2 ;
                    p1 = repmat(elems(:,1),[nDiv 1]) ;
                    p2 = elems(:,2:end-1) ;
                    p3 = elems(:,3:end) ;
                    elems = [p1(:) p2(:) p3(:)] ;
                end
            % Triangulation class constructor
                elems = triangulation(elems,refPoints) ; 
            % Initialization
                indPts = ones(size(Points,1),3) ;
            % Enclosing element and barycentric positions
                [elmt,weights] = pointLocation(elems,Points) ; 
            % Outside-of-old-mesh newPoints need to be fixed (elmt=NaN)
                outside = isnan(elmt) ;
                if any(outside)
                    switch extrap
                        case 'extrap'
                            % Outside Points
                                outsidePts = Points(outside,:) ;
                            % Find the closest element
                                C = circumcenter(tri) ;
                                [~,clstElmt] = min(sum((reshape(outsidePts,[],1,2)-reshape(C,1,[],2)).^2,3),[],2) ;
                                elmt(outside) = clstElmt ;
                            % Find the (extended) coordinates of the new point in each old closest element frame
                                elmtNodes = elems(clstElmt,:) ;
                                p1 = refPoints(elmtNodes(:,1),:)' ; 
                                p2 = refPoints(elmtNodes(:,2),:)' ; 
                                p3 = refPoints(elmtNodes(:,3),:)' ; 
                                v1 = p2-p1 ; v2 = p3-p1 ; % element's frame vectors
                                m = [] ; % coordinates in this frame
                                for pp = 1:length(clstElmt)
                                    m(pp,:) = [v1(:,pp) v2(:,pp)]\(outsidePts(pp,:)'-p1(:,pp)) ;
                                end
                            % Assign the weights
                                indPts(outside,:) = elmtNodes ;
                                weights(outside,:) = [1-m(:,1)-m(:,2) m(:,1) m(:,2)] ;
                        otherwise
                            weights(outside,:) = extrap ;
                    end
                end
            % Inside nodes
                indPts(~outside,:) = elems(elmt(~outside),:) ;
            % Transfer Matrix
                nodes = (1:size(Points,1))'*[1 1 1] ;
                T = sparse(nodes(:),indPts(:),weights(:),size(Points,1),size(refPoints,1)) ;
        end
    end


% ------- DATA DIFFERENTIATION AND INTEGRATION --------------------------

    function [D1,D2] = gradMat(obj,refPoints,onNodes)
    % Return the differentiation matrices so that df_dxi = Di*f(:)
    % size(Di) = [nElems nNodes] : f is defined on nodes, the gradient
    % is constant over elements
        if nargin==1 ; refPoints = obj.Points ; end
        if nargin<3 ; onNodes = obj.DataOnNodes ; end
        if size(obj.Elems,2)==3 % TRIANGLES ONLY
            nTris = size(obj.Triangles,1) ; 
            nPts = size(refPoints,1) ;
            nDims = size(refPoints,2) ;
            nNodesByElem = size(obj.Triangles,2) ; %(==3) here
            % Face points 
                Polygons = reshape(refPoints(obj.Triangles(:),:),nTris,[],nDims) ;
            % Areas
                Areas = polyarea(Polygons(:,:,1)',Polygons(:,:,2)').' ;
            % Derivatives
                % Triangles shape indicators
                    %aa = cross(Polygons(:,:,1),Polygons(:,:,2),2) ;
                    bb = circshift(Polygons(:,:,2),2,2) - circshift(Polygons(:,:,2),1,2) ;
                    cc = circshift(Polygons(:,:,1),1,2) - circshift(Polygons(:,:,1),2,2) ;
                % Sparse matrices indices
                    iii = repmat((1:nTris)',[1 nNodesByElem]) ;
                    jjj = obj.Triangles ;
                    vvvD1 = .5.*bb./Areas(:) ;
                    vvvD2 = .5.*cc./Areas(:) ;
                % Build sparse matrices
                    D1 = sparse(iii(:),jjj(:),vvvD1(:),nTris,nPts) ;
                    D2 = sparse(iii(:),jjj(:),vvvD2(:),nTris,nPts) ;
        else % THERE IS DIFFERENT ELEMENT TYPES
            % Build the element matrices
                elems = obj.Elems ;
            % Coordinates
                X = refPoints(:,1) ; Y = refPoints(:,2) ;
            % Compute the mean gradient
                eee = [] ;
                nnn = [] ;
                vvvD1 = [] ;
                vvvD2 = [] ;
                for elmt = 1:size(elems,1)
                    % Element topology
                        nodes = elems(elmt,~isnan(elems(elmt,:)))' ;
                        area = polyarea(X(nodes),Y(nodes)) ;
                        edges = [nodes circshift(nodes,1,1)] ;
                    % Normals
                        normals = [-diff(Y(edges),1,2) diff(X(edges),1,2)] ;
                        dl = sum(normals.^2,2) ; 
                        midPts = [mean(X(edges),2) mean(Y(edges),2)] ;
                        switched = inpolygon(midPts(:,1)+normals(:,1)*0.01,midPts(:,2)+normals(:,2)*0.01,X(nodes),Y(nodes)) ;
                        normals(switched,:) = -normals(switched,:) ;
                        normals = normals./dl ;
                    % Sparse matrices indices
                        eee = [eee ; elmt*ones(numel(edges),1)] ;
                        nnn = [nnn ; edges(:)] ;
                        vvvD1 = [vvvD1 ; reshape(1/2/area.*dl.*normals(:,1).*[1 1],[],1)] ;
                        vvvD2 = [vvvD2 ; reshape(1/2/area.*dl.*normals(:,2).*[1 1],[],1)] ;
                end
            % Differentiation matrices
                D1 = sparse(eee(:),nnn(:),vvvD1(:),size(elems,1),size(refPoints,1)) ;
                D2 = sparse(eee(:),nnn(:),vvvD2(:),size(elems,1),size(refPoints,1)) ;
        end
    % Apply to nodes if needed (mean over connected elements)
        if onNodes
            T = obj.elem2nod ;
            D1 = T*D1 ; D2 = T*D2 ;
        end
    end

    function [v,areas] = integMat(obj,refPoints)
    % Return a vector v so that \int(f(x).dx) = v*f
    % f can be defined on nodes (DataOnNodes = true) or on elements
    % (DataOnNodes = false)
        if nargin==1 ; refPoints = obj.Points ; end
        elems = obj.Elems ;
        % Element Areas
            areas = zeros(size(obj.Elems,1),1) ;
            for nNodes = 1:size(elems,2) % For each polygonal element shape..
                hasNNodes = sum(~isnan(elems),2)==nNodes ;
                if ~any(hasNNodes) ; continue ; end
                x = reshape(refPoints(elems(hasNNodes,1:nNodes),:),[],nNodes,2) ;
                areas(hasNNodes) = polyarea(x(:,:,1),x(:,:,2),2) ;
            end
        % Integration matrix
            if obj.DataOnNodes
                elem2nod = obj.elem2nod('ones') ;
                elem2nod = elem2nod./sum(elem2nod,1) ;
                v = areas(:)'*elem2nod' ;
            else
                v = areas(:)' ;
            end
    end
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH DATA HANDLING AND COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    function set.DataOnNodes(obj,val)
    % Put every data on nodes
        obj.DataOnNodes = val ;
        obj.computeDataFields ;
    end

    function DATA = computeDataFields(obj,onNodes,fr)
    % Compute all (scalar) data fields associated to the object motion
    % Empty fields are here as section headings
        if isempty(obj.MovingPoints) ; return ; end
        if nargin<2 || isempty(onNodes) ; onNodes = obj.DataOnNodes ; end
        [nPoints,nCoord,nFrames] = size(obj.MovingPoints) ;
        if nargin<3 ; fr = 1:nFrames ; end
        % Init the structure
            computeAll = isempty(obj.DataFields) || nargin<3 ;
            DATA = struct() ;
        % REFERENCE CONFIGURATION
            if computeAll % re-determine the reference config
                DATA.Reference = 'Reference' ;
                X = reshape(permute(obj.MovingPoints,[1 3 2]),[nPoints*nFrames nCoord]) ;
                % Indice of the first valid frame (for mechanical fields)
                    [~,firstValidFrame] = max(all(~isnan(obj.MovingPoints),2),[],3) ;
                    DATA.FirstValidFrame = firstValidFrame(:) ;
                    indFirstValid = sub2ind([nPoints nFrames],1:nPoints,firstValidFrame(:)') ;
                    DATA.X1first = X(indFirstValid,1) ;
                    DATA.X2first = X(indFirstValid,2) ;
                % Indice of the reference frame for gradient data computation
                % taken as the first frame for which all the points are valid
                    dataRef = max(DATA.FirstValidFrame).*ones(nPoints,1) ;
                    DATA.StrainReference = dataRef ;
                    indRef = sub2ind([nPoints nFrames],1:nPoints,dataRef') ;
                    DATA.X1ref = X(indRef,1) ;
                    DATA.X2ref = X(indRef,2) ;
            else % copy from the existing data
                toCopy = fieldnames(obj.DataFields)' ;
                toCopy = toCopy(1:7) ;
                for field = string(toCopy)
                    DATA.(field) = obj.DataFields.(field{1}) ;
                end
            end
        % POSITION
            DATA.Position = 'Position' ;
            % NaNs
                DATA.NaN = obj.MovingPoints(:,1,fr)*NaN ;
            % Current Coordinates
                DATA.x1 = obj.MovingPoints(:,1,fr) ;
                DATA.x2 = obj.MovingPoints(:,2,fr) ;
        % Displacements
            DATA.Displacement = 'Displacement' ;
            DATA.u1 = DATA.x1 - DATA.X1first ;
            DATA.u2 = DATA.x2 - DATA.X2first ;
            DATA.U = sqrt(DATA.u1.^2 + DATA.u2.^2) ;
        % Velocity (pixels/image)
            DATA.Velocity = 'Velocity' ; 
            DATA.v1 = cat(3,zeros(nPoints,1),diff(DATA.u1,1,3)) ;
            DATA.v2 = cat(3,zeros(nPoints,1),diff(DATA.u2,1,3)) ;
            DATA.V = sqrt(DATA.v1.^2 + DATA.v2.^2) ;
        % DATA FIELDS WITH GRADIENT
            if ~isempty(obj.Triangles) || ~isempty(obj.Quadrangles)
                DATA = surfaceDataFields(obj,DATA,onNodes) ;
            elseif size(obj.Elems,2)==2
                DATA = barDataFields(obj,DATA) ;
            end
        % SAVE IN THE OBJECT
            if computeAll
                obj.DataFields = DATA ;
            else
                toCopy = setdiff(fieldnames(obj.DataFields)',toCopy) ;
                for field = string(toCopy)
                    if ischar(obj.DataFields.(field)) ; continue ; end
                    obj.DataFields.(field)(:,:,fr) = DATA.(field) ;
                end
            end
    end
    
    function DATA = surfaceDataFields(obj,DATA,onNodes)
    % Compute data fields associated to surfacic elements
    % Derivation matrices
        [nPoints,~,nFrames] = size(DATA.x1) ;
        [D1,D2] = gradMat(obj,[DATA.X1ref DATA.X2ref],onNodes) ;
    % Transformation Gradient
        DATA.Transformation = 'Transformation' ; 
        x1 = DATA.x1 ;
        x2 = DATA.x2 ;
        DATA.F11 = reshape(D1*squeeze(x1),[],1,nFrames) ;
        DATA.F12 = reshape(D2*squeeze(x1),[],1,nFrames) ;
        DATA.F21 = reshape(D1*squeeze(x2),[],1,nFrames) ;
        DATA.F22 = reshape(D2*squeeze(x2),[],1,nFrames) ;
    % Modify the transformation if needed (reference change)
    % F(t) = dx(t)/dXfirst 
    %      = dx(t)/dXref * dXref/dXfirst
    %      = Fref(t)*(Fref(tfirst))^{-1}
        changeRef = DATA.FirstValidFrame~=DATA.StrainReference ;
        if any(changeRef)
            DATA.Fref11 = DATA.F11 ;
            DATA.Fref12 = DATA.F12 ;
            DATA.Fref21 = DATA.F21 ;
            DATA.Fref22 = DATA.F22 ;
        % First valid linear index
            nData = size(DATA.F11,1) ;
            [~,fvi] = max(~isnan(DATA.F11),[],3) ;
            fvi = sub2ind([nData nFrames],(1:nData)',fvi) ;
        % Inverse transformation (Fref(tfirst))^{-1}
            [iF11,iF22,iF21,iF12] = inverse(obj,DATA.F11(fvi),DATA.F22(fvi),DATA.F21(fvi),DATA.F12(fvi)) ;
        % New Values
            DATA.F11 = DATA.Fref11.*iF11 + DATA.Fref12.*iF21 ;
            DATA.F12 = DATA.Fref11.*iF12 + DATA.Fref12.*iF22 ;
            DATA.F21 = DATA.Fref21.*iF11 + DATA.Fref22.*iF21 ;
            DATA.F22 = DATA.Fref21.*iF12 + DATA.Fref22.*iF22 ;
        end
    % Other transfoermation data
        DATA.J2D = DATA.F11.*DATA.F22 - DATA.F12.*DATA.F21 ; % 2D Jacobian
        DATA.J = DATA.J2D*0+1 ; % 3D Jacobian=1 (incompressible)
        DATA.F33 = 1./(DATA.J2D) ; %J = det(F_3D) = det(F_2D)*F33 = 1
    % Cauchy Strains
        DATA.CauchyStrains = 'Cauchy Strains' ; 
        DATA.C11 = DATA.F11.^2 + DATA.F21.^2 ;
        DATA.C22 = DATA.F12.^2 + DATA.F22.^2 ;
        DATA.C12 = DATA.F11.*DATA.F12 + DATA.F21.*DATA.F22 ;
        DATA.C33 = DATA.F33.^2 ;
        [l1,l2,~,DATA.Ctheta] = eigenValues(obj,DATA.C11,DATA.C22,DATA.C12) ;
        DATA.lambda1 = sqrt(l1) ; 
        DATA.lambda2 = sqrt(l2) ; 
        DATA.lambda3 = 1./(DATA.lambda1.*DATA.lambda2) ;
    % Green-Lagrange Strains
        DATA.Green_LagrangeStrains = 'Green-Lagrange Strains' ; 
        DATA.L11 = 0.5*(DATA.C11 - 1) ;
        DATA.L22 = 0.5*(DATA.C22 - 1) ;
        DATA.L12 = 0.5*(DATA.C12) ;
        DATA.L33 = 0.5*(DATA.C33 - 1) ;
        [DATA.Ldev11,DATA.Ldev22,DATA.Ldev33,DATA.Ldev12,DATA.Lmean,DATA.Leq] = deviatoric(obj,DATA.L11,DATA.L22,DATA.L33,DATA.L12) ;
        [DATA.Le1,DATA.Le2,DATA.Ltau,DATA.Ltheta] = eigenValues(obj,DATA.L11,DATA.L22,DATA.L12) ;
    % Strain rate tensor (using Lagrangian fields)
    % D = inv(F').Ldot.inv(F)
        Ldot11 = cat(3,0*DATA.L11(:,:,1),diff(DATA.L11,1,3)) ;
        Ldot22 = cat(3,0*DATA.L22(:,:,1),diff(DATA.L22,1,3)) ;
        Ldot12 = cat(3,0*DATA.L12(:,:,1),diff(DATA.L12,1,3)) ;
        idet2 = (DATA.F11.*DATA.F22 - DATA.F21.*DATA.F12).^-2 ;
        DATA.StrainRate = 'Strain Rate' ; 
        DATA.D11 = idet2.*(Ldot11.*DATA.F22.^2 + Ldot22.*DATA.F21.^2 - 2*Ldot12.*DATA.F21.*DATA.F22) ;
        DATA.D22 = idet2.*(Ldot11.*DATA.F12.^2 + Ldot22.*DATA.F11.^2 - 2*Ldot12.*DATA.F12.*DATA.F11) ;
        DATA.D12 = idet2.*(Ldot12.*(DATA.F11.*DATA.F22 + DATA.F21.*DATA.F12) - Ldot11.*DATA.F22.*DATA.F12 - Ldot22.*DATA.F11.*DATA.F21) ;
        DATA.D33 = -DATA.D11-DATA.D22 ;
        [DATA.Ddev11,DATA.Ddev22,DATA.Ddev33,DATA.Ddev12,DATA.Dmean,DATA.Deq] = deviatoric(obj,DATA.D11,DATA.D22,DATA.D33,DATA.D12) ;
        [DATA.De1,DATA.De2,DATA.Dtau,DATA.Dtheta] = eigenValues(obj,DATA.D11,DATA.D22,DATA.D12) ;
    % True Strains
        DATA.TrueStrain = 'True Strain' ; 
        DATA.TS11 = cumsum(DATA.D11,3) ;
        DATA.TS22 = cumsum(DATA.D22,3) ;
        DATA.TS12 = cumsum(DATA.D12,3) ;
        DATA.TS33 = cumsum(DATA.D33,3) ;
        [DATA.TSdev11,DATA.TSdev22,DATA.TSdev33,DATA.TSdev12,DATA.TSmean,DATA.TSeq] = deviatoric(obj,DATA.TS11,DATA.TS22,DATA.TS33,DATA.TS12) ;
        [DATA.TSe1,DATA.TSe2,DATA.TStau,DATA.TStheta] = eigenValues(obj,DATA.TS11,DATA.TS22,DATA.TS12) ;
    % Linearized Green-Lagrange Strains
        DATA.LinearizedStrains = 'Linearized Strains' ; 
        DATA.E11 = reshape(D1*squeeze(DATA.u1),[],1,nFrames) ;
        DATA.E22 = reshape(D2*squeeze(DATA.u2),[],1,nFrames) ;
        DATA.E12 = 0.5*(reshape(D2*squeeze(DATA.u1),[],1,nFrames) + reshape(D1*squeeze(DATA.u2),[],1,nFrames)) ;
        [DATA.Ee1,DATA.Ee2,DATA.Etau,DATA.Etheta] = eigenValues(obj,DATA.E11,DATA.E22,DATA.E12) ;
    end
    
    function DATA = barDataFields(obj,DATA)
    % Compute data fields associated to rod elements
        nPoints = size(obj.Points,1) ;
        bars = obj.Elems(:,1:2) ;
        nBars = size(bars,1) ;
        nFrames = size(obj.MovingPoints,3) ;
    % Differenciation matrix
        D = sparse((1:nBars)'*[1 1],bars,ones(nBars,1).*[-1 1],nBars,nPoints) ;
    % Differetiation position
        X = reshape(cat(2,DATA.x1,DATA.x2),nPoints,2*nFrames) ;
        DX = reshape(D*X,[],2,nFrames) ;
    % Bars as a complex number
        P = DX(:,1,:)+1i*DX(:,2,:) ;
    % Length Measurement
        DATA.Length = 'Length' ; 
        DATA.L = abs(P) ;
        DATA.dL = cat(3,zeros(size(DX,1),1),diff(DATA.L,1,3)) ;
        DATA.dLtot = cumsum(DATA.dL,3) ;
    % Strains
        DATA.Length = 'Stretch' ; 
        DATA.E = DATA.dLtot./DATA.L(:,1) ;
        DATA.lambda = DATA.E + 1 ;
        DATA.D = DATA.dL./DATA.L ;
        DATA.TS = cumsum(DATA.D,3) ;
    % Edge Rotation
        DATA.Rotations = 'Rotation' ; 
        DATA.A = angle(P) ;
        DATA.dA = angle(cat(3,zeros(size(DX,1),1),P(:,:,2:end)./P(:,:,1:end-1))) ;
        DATA.dAtot = cumsum(DATA.dA,3) ;
    end

    function [e1,e2,tau,theta] = eigenValues(obj,M11,M22,M12)
    % Return eigen values of a 2x2 symmetric matrix
        aaa = 1/2*(M11+M22) ;
        bbb = sqrt(1/4*(M11-M22).^2 + M12.^2) ;
        e1 = aaa+bbb ; % Major
        e2 = aaa-bbb ; % Minor
        tau = bbb ; % Shear
        theta = 1/2*atan(2*M12./(M11-M22)) ; % Principal Angle
    end

    function [iM11,iM22,iM21,iM12] = inverse(obj,M11,M22,M21,M12)
    % Return the inverse of a 2x2 matrix
        if nargin<5 ; M12 = M21 ; end
        det = M11.*M22-M21.*M12 ;
        iM11 = M22./det ;
        iM22 = M11./det ;
        iM21 = -M21./det ;
        iM12 = -M12./det ;
    end

    function [Mdev11,Mdev22,Mdev33,Mdev12,Mmean,Meq] = deviatoric(obj,M11,M22,M33,M12)
    % Return the deviatoric part of a symmetric tensor
        Mmean = 1/3*(M11+M22+M33) ;
        Mdev11 = M11-Mmean ;
        Mdev22 = M22-Mmean ;
        Mdev33 = M33-Mmean ;
        Mdev12 = M12 ;
        Meq = sqrt(2/3*(Mdev11.^2 + Mdev22.^2 + Mdev33.^2 + 2*Mdev12.^2)) ;
    end
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEED PREVIEW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    function initSeedPreview(obj,ax)
        % INIT THE MESH
            axes(ax) ;
            triMesh = patch(obj.Points(:,1),obj.Points(:,2),NaN*obj.Points(:,2) ...
                            ,'vertices',obj.Points...
                            ,'faces',obj.Elems...
                            ,'edgecolor','k'...
                            ,'EdgeAlpha', 0.2 ...
                            ,'linewidth',0.1...
                            ,'facecolor','interp'...
                            ,'facealpha',0.5 ...
                            ,'hittest','off' ...
                            ,'Visible','off' ...
                            ,'tag',obj.Name ...
                            ,'DisplayName','mesh' ...
                            ) ;
            scatt = scatter(obj.Points(:,1),obj.Points(:,2)...
                            ,200 ... scatter size
                            ,'hittest','off' ...
                            ,'marker','.' ...
                            ,'Visible','off' ...
                            ,'tag',obj.Name ...
                            ,'DisplayName','scatter' ...
                            ) ;
            contour = patch(obj.Points(:,1),obj.Points(:,2),NaN*obj.Points(:,2) ...
                            ,'edgecolor','interp'...
                            ,'EdgeAlpha', 1 ...
                            ,'linewidth',1.5 ...
                            ,'hittest','off' ...
                            ,'Visible','off' ...
                            ,'tag',obj.Name ...
                            ,'DisplayName','contour' ...
                            ) ;
            edges = patch(obj.Points(:,1),obj.Points(:,2),NaN*obj.Points(:,2) ...
                            ,'vertices',obj.Points...
                            ,'faces',obj.Edges...
                            ,'edgecolor','flat'...
                            ,'EdgeAlpha', 1 ...
                            ,'linewidth',2 ...
                            ,'facecolor','none'...
                            ,'hittest','off' ...
                            ,'Visible','off' ...
                            ,'tag',obj.Name ...
                            ,'DisplayName','edges' ...
                            ) ;
        % ADD A COLORBAR
            clrbr = colorbar(ax) ;
            fig = ax.Parent ;
            switch fig.Position(3)<fig.Position(4)
                case true
                    clrbr.Location = 'east' ;
                case false
                    clrbr.Location = 'north' ;
            end
            clrbr.Ruler.Exponent = 0 ;
            clrbr.FontWeight = 'bold' ;
            if mean(obj.refImgs{1}(:))/max(getrangefromclass(obj.refImgs{1}))<0.5
                clrbr.Color = 'w' ;
            end
        % ADD THE MENU BAR
            submenus = gobjects(0) ;
            % Data to plot
                mData = uimenu(ax.Parent,'Label','Data') ;
                % Re-compute DataFields if needed
                    if isempty(obj.DataFields) || isempty(fieldnames(obj.DataFields)) ; obj.computeDataFields() ; end
                % Fill Menus
                    if ~isempty(fieldnames(obj.DataFields))
                        % Get dete fields
                            fieldNames = fieldnames(obj.DataFields) ;
                        % Add a menu for each of them
                            letters = @(str)str(uint8(str)>=65 & uint8(str)<=122) ;
                            for f = 1:numel(fieldNames)
                                % Add a menu if needed (the field is empty)
                                    if ischar(obj.DataFields.(fieldNames{f}))
                                        mParent = uimenu(mData,'Label',obj.DataFields.(fieldNames{f})) ;
                                        continue ;
                                    end
                                % Add the submenu
                                    submenus(end+1) = uimenu(mParent,'Label',fieldNames{f}) ;
                                % Add a separator if needed (one letter different)
                                    if numel(mParent.Children)>1 && ~strcmpi(letters(fieldNames{f}),letters(fieldNames{f-1}))
                                        submenus(end).Separator = 'on' ;
                                    end
                            end
                        % Check the default choice
                            submenus(ismember({submenus.Label},'U')).Checked = 'on' ;
                    end
            % DISPLAY
                mDisplay = uimenu(ax.Parent,'Label','Display') ;
                % Plot Type
                    mPlot = uimenu(mDisplay,'Label','Plot') ;
                        submenus(end+1) = uimenu(mPlot,'Label','Mesh','checked','on') ;
                        submenus(end+1) = uimenu(mPlot,'Label','Contours') ;
                        submenus(end+1) = uimenu(mPlot,'Label','Scatter') ;
                        submenus(end+1) = uimenu(mPlot,'Label','Edges') ;
                % Data Scale
                    mScale = uimenu(mDisplay,'Label','Scale') ;
                        submenus(end+1) = uimenu(mScale,'Label','Linear','checked','on') ;
                        submenus(end+1) = uimenu(mScale,'Label','Log10') ;
                        submenus(end+1) = uimenu(mScale,'Label','Cummul') ;
                % Color Scale
                    mColors = uimenu(mDisplay,'Label','Colors') ;
                        mClrFrames = uimenu(mColors,'Label','Frames') ;
                            submenus(end+1) = uimenu(mClrFrames,'Label','Current','checked','on') ;
                            submenus(end+1) = uimenu(mClrFrames,'Label','All') ;
                        mClrLims = uimenu(mColors,'Label','Limits') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','min-max','checked','on') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','0-max') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','symmetric') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','1*sigma','separator','on') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','2*sigma') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','3*sigma') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','4*sigma') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','5*sigma') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','6*sigma') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','cummul','separator','on') ;
                        mClrSteps = uimenu(mColors,'Label','Steps') ;
                            for ss = 4:4:24 ; submenus(end+1) = uimenu(mClrSteps,'Label',num2str(ss)) ; end
                            submenus(end+1) = uimenu(mClrSteps,'Label','Continuous','checked','on') ;
                        submenus(end+1) = uimenu(mColors,'Label','Custom') ;
                % Axes Behavior
                    mAxesMode = uimenu(mDisplay,'Label','Axes') ;
                        submenus(end+1) = uimenu(mAxesMode,'Label','Fixed','checked','on') ;
                        submenus(end+1) = uimenu(mAxesMode,'Label','Follow') ;
            % Common Properties
                set(submenus,'callback',@(src,evt)obj.updateSeedMenus(src,ax)) ;
            % UserData in axes to choose the data to plot
                ax.UserData.dataLabel = 'U' ;
                ax.UserData.plotType = 'Mesh' ;
                ax.UserData.dataScale = 'Linear' ;
                ax.UserData.clrMode = 'Preset' ;
                ax.UserData.clrFramesLabel = 'Current' ;
                ax.UserData.clrLimsLabel = 'min-max' ;
                ax.UserData.clrStepsLabel = 'Continuous' ;
                ax.UserData.axesPositionReference = [] ; % [frame posX posY]
    end


    function updateSeedMenus(obj,subMenu,ax)
        % Uncheck all subMenuItems
            set(subMenu.Parent.Children,'checked','off')
        % Check the selected item
            subMenu.Checked = 'on' ;
        % Change the ax userData if needed
            switch subMenu.Parent.Label
                case 'Frames'
                    ax.UserData.clrFramesLabel = subMenu.Label ;
                    ax.UserData.clrMode = 'Preset' ;
                case 'Limits'
                    ax.UserData.clrLimsLabel = subMenu.Label ;
                    ax.UserData.clrMode = 'Preset' ;
                case 'Plot'
                    ax.UserData.plotType = subMenu.Label ;
                case 'Scale'
                    ax.UserData.dataScale = subMenu.Label ;
                case 'Steps'
                    ax.UserData.clrStepsLabel = subMenu.Label ;
                    ax.UserData.clrMode = 'Preset' ;
                case 'Axes'
                    switch subMenu.Label
                        case 'Fixed'
                            ax.UserData.axesPositionReference = [] ;
                        case 'Follow'
                            % The axes position will start following the current position of the axes center
                                cen = [mean(ax.XLim) mean(ax.YLim)] ; ran = [range(ax.XLim) range(ax.YLim)] ;
                                frame = ax.UserData.currentFrame ;
                            % Position of the current point in the reference config.
                                cen = full(obj.interpMat(cen,'extrap',obj.MovingPoints(:,:,frame))*obj.Points) ;
                            % Record
                                ax.UserData.axesPositionReference = [cen ran] ; % [frame posX posY]
                                disp([obj.Name ' preview will follow with the position reference ' mat2str(ax.UserData.axesPositionReference)]) ;
                    end
                case 'Colors' 
                    switch subMenu.Label
                        case 'Custom'
                            % Let the user choose the color scale
                                prompt = {'Min:','Max:','Steps:'};
                                dlgtitle = 'Custom Color Scale Parameters';
                                dims = [1 35];
                                definput = {num2str(min(caxis(ax))),num2str(max(caxis(ax))),num2str(size(colormap(ax),1))};
                                answer = inputdlg(prompt,dlgtitle,dims,definput) ;
                                if isempty(answer) ; return ; end
                            % Set the axis color properties
                                caxis(ax,[str2num(answer{1}),str2num(answer{2})])
                                colormap(ax,jet(str2num(answer{3}))) ;
                            % Turn on custom mode
                                ax.UserData.clrMode = 'Custom' ;
                    end
                otherwise % Data Menus
                    ax.UserData.dataLabel = subMenu.Label ;
                    % Uncheck all Data Menus
                        set([subMenu.Parent.Parent.Children],'checked','off')
                        set(cat(1,subMenu.Parent.Parent.Children.Children),'checked','off')
                    % Check the parent menu
                        subMenu.Parent.Checked = 'on' ;
                        subMenu.Checked = 'on' ;
            end
        % Update the preview
            updateSeedPreview(obj,[],ax)
    end


    function updateSeedPreview(obj,~,ax)
        % Has the preview been initialized ?
            % Search for the object
                gH = findobj(ax,'tag',obj.Name) ;
            % If it is not found, re-init
                if isempty(gH)
                    initSeedPreview(obj,ax) ;
                    gH = findobj(ax,'tag',obj.Name) ;
                end
        % Structure containing plot Data
            PlotStruct = [] ;
        % Get the frame
            CurrentFrame = ax.UserData.currentFrame ;
        % If the frame is not valid, return
            if ~(CurrentFrame>0 && CurrentFrame<=size(obj.MovingPoints,3)) ; return ; end
        % Current configuration
            PlotStruct.Vertices = obj.MovingPoints(:,:,CurrentFrame) ;
        % Data Field to plot
            Data = [] ;
            % Re-compute data fields if needed
                if ~isfield(obj.DataFields,ax.UserData.dataLabel) ; obj.computeDataFields ; end
            % Extract the data
                Data = obj.DataFields.(ax.UserData.dataLabel) ;
            % If no data, return
                if isempty(Data) ; return ; end
            % Re-format
                Data = permute(Data,[1 3 2]) ;
        % Data Transformation
            switch ax.UserData.dataScale
                case 'Linear' % Do nothing
                    ax.ColorScale = 'linear' ;
                case 'Log10'
                    ax.ColorScale = 'log' ;
                    %Data = log10(abs(Data)) ;
                case 'Cummul'
                    [sortData,sortInd] = sort(Data,1) ;
                    sortInd = sortInd + ((1:size(Data,2))-1)*size(Data,1) ;
                    newData = cumsum(sortData,1) ;
                    newData = newData./newData(end,:) ;
                    newIndices(sortInd(:)) = 1:numel(Data) ;
                    Data = reshape(newData(newIndices),size(Data)) ;
            end
        % Extract the Current Data
            CurrentFrame = min(size(Data,2),CurrentFrame) ; % Allow a field to be constant
            PlotStruct.CData = Data(:,CurrentFrame) ;
        % COLORS
            if strcmp(ax.UserData.clrMode,'Preset')
                % Color Frames
                    switch ax.UserData.clrFramesLabel
                        case 'Current'
                            colorData = Data(:,CurrentFrame) ;
                        case 'All'
                            colorData = Data ;
                    end
                    minData = min(colorData(:)) ;
                    maxData = max(colorData(:)) ;
                % Color Limits
                    CLim = [minData maxData] ;
                    if ~any(isnan(CLim))
                        % Set CAXIS
                            switch ax.UserData.clrLimsLabel
                                case '0-max'
                                    if sign(minData)==sign(maxData)
                                        if sign(minData)==-1
                                            CLim = [minData 0] ;
                                        else
                                            CLim = [0 maxData] ;
                                        end
                                    else
                                        CLim = [minData maxData] ;
                                    end
                                case 'min-max'
                                    CLim = [minData maxData] ;
                                case 'symmetric'
                                    CLim = max(abs([minData maxData]))*[-1 1] ;
                                otherwise % N*sigma
                                    % Get the factor N
                                        N = str2double(ax.UserData.clrLimsLabel(1)) ;
                                    % Compute mean and variances
                                        valid = ~isnan(colorData(:)) ;
                                        avg = mean(colorData(valid)) ;
                                        ec = std(colorData(valid)) ;
                                    % Set Color Limits
                                        CLim = min(max(avg+N*ec*[-1 1],minData),maxData) ;
                            end
                            if range(CLim)<eps ; CLim = CLim(1)+abs(CLim(1))*[-1 1]*eps ; end
                            if any(CLim~=caxis(ax)) ; caxis(ax,CLim) ; end
                    end
                % Color Steps
                    steps = ax.UserData.clrStepsLabel ;
                    switch steps
                        case 'Continuous'
                            steps = 1000 ;
                        otherwise
                            steps = str2num(steps) ;
                    end
                    if steps~=size(colormap(ax),1) ; colormap(ax,jet(steps)) ; end
            else
                CLim = caxis(ax) ;
                steps = size(colormap(ax),1) ;
            end
        % UPDATE THE DISPLAY OBJECT
            % Set objects not visible
                set(gH,'Visible','off') ;
            % Depending on the kind of plot..
                switch ax.UserData.plotType
                    case 'Mesh'
                        gHi = findobj(gH,'DisplayName','mesh') ;
                        % Apply data to mesh
                            if size(PlotStruct.CData,1) == size(obj.Elems,1)
                                gHi.FaceColor = 'flat' ;
                            else
                                gHi.FaceColor = 'interp' ;
                            end
                        % Display
                            gHi.Vertices = PlotStruct.Vertices ;
                            gHi.FaceVertexCData = PlotStruct.CData ;
                    case 'Contours'
                        gHi = findobj(gH,'DisplayName','contour') ;
                        % Force data on nodes
                            if ~obj.DataOnNodes && size(PlotStruct.CData,1)==size(obj.Elems,1)
                                PlotStruct.CData = obj.elem2nod*PlotStruct.CData ; 
                            end
                        % Contour levels
                            nL = min(40,steps) ;
                            lvl = linspace(CLim(1),CLim(2),2*nL+1) ;
                            lvl = lvl(2:2:end) ;
                            lvl(isnan(lvl)) = [] ;
                        % Compute contour lines
                            C = NaN*[1 1 1] ;
                            if ~isempty(lvl)
                                C = tricontours(obj.Triangles,PlotStruct.Vertices,PlotStruct.CData,lvl) ;
                                C = reshape(C,[],3) ; % A list of vertices with interleaved NaNs
                            end
                        % Display
                            gHi.Vertices = C.*[1 1 0] ;
                            gHi.CData = C(:,3) ;
                            gHi.Faces = 1:size(gHi.Vertices,1) ;
                    case 'Scatter'
                        gHi = findobj(gH,'DisplayName','scatter') ;
                        % Force data on nodes
                            if obj.DataOnNodes && size(PlotStruct.CData,1)==size(obj.Elems,1)
                                PlotStruct.CData = obj.elem2nod*PlotStruct.CData ; 
                            end
                        % Display
                            gHi.XData = PlotStruct.Vertices(:,1) ;
                            gHi.YData = PlotStruct.Vertices(:,2) ;
                            gHi.CData = PlotStruct.CData ;
                    case 'Edges'
                        gHi = findobj(gH,'DisplayName','edges') ;
                        % Force data on nodes
                            if size(PlotStruct.CData,1)==size(obj.Elems,1)
                                PlotStruct.CData = obj.elem2nod*PlotStruct.CData ; 
                            end
                        % Display
                            gHi.Vertices = PlotStruct.Vertices ;
                            gHi.FaceVertexCData = PlotStruct.CData ;
                end
            % Make the object visible
                gHi.Visible = 'on' ;
        % AXES LIMITS
            % Follow point if needed
                if ~isempty(ax.UserData.axesPositionReference) && ~any(isnan(ax.UserData.axesPositionReference))
                    % Center and range of axis limits
                        cen = (obj.interpMat(ax.UserData.axesPositionReference(1:2))*obj.MovingPoints(:,:,ax.UserData.currentFrame)) ; 
                        ran = (ax.UserData.axesPositionReference(3:4)) ;
                    % Update limits
                        if ~any(isnan([cen ran]))
                            [nI,nJ] = size(obj.refImgs{1}) ;
                            ax.XLim = max(0.5,min(nJ+0.5,cen(1)+ran(1)*[-1 1]/2)) ;
                            ax.YLim = max(0.5,min(nI+0.5,cen(2)+ran(2)*[-1 1]/2)) ;
                        end
                end
    end
        
        
end
    
end 