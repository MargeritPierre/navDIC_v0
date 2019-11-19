classdef navDICSeed_2D_DistMesh < navDICSeed_2D_Surface
   
    properties
        Triangles = [] ;
        Quadrangles = [] ;
        Edges = [] ;
        ROI = [] ;
        Shapes = [] ;
        DataOnNodes = false ;
    end
    
    properties (Transient)
        DataFields = struct() ;
    end
    
    methods
        
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
%                 obj.Displacements = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
%                 obj.Strains = ones(size(obj.Points,1),3,hd.nFrames)*NaN ;
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
        
        
        function T = interpMat(obj,Points,extrap,refPoints)
        % Return the interpolation matrix of the mesh on query points so
        % that U(pts) = T*U(Nodes)
            if nargin<3 ; extrap = 'extrap' ; end
            if nargin<4 ; refPoints = obj.Points ; end
            % Triangulation class constructor
                tri = triangulation(obj.Triangles,refPoints) ; 
            % Initialization
                indPts = ones(size(Points,1),3) ;
            % Enclosing element and barycentric positions
                [elmt,weights] = pointLocation(tri,Points) ; 
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
                                elmtNodes = obj.Triangles(clstElmt,:) ;
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
                indPts(~outside,:) = obj.Triangles(elmt(~outside),:) ;
            % Transfer Matrix
                nodes = (1:size(Points,1))'*[1 1 1] ;
                T = sparse(nodes(:),indPts(:),weights(:),size(Points,1),size(refPoints,1)) ;
        end
        
        
        function [D1,D2] = gradMat(obj,refPoints)
        % Return the differentiation matrices so that df_dxi = Di*f(:)
        % size(Di) = [nElems nNodes] : f is defined on nodes, the gradient
        % is constant over elements
            if nargin==1 ; refPoints = obj.Points ; end
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
            % Apply to nodes if needed (mean over connected triangles)
                if obj.DataOnNodes
                    T = obj.tri2nod ;
                    D1 = T*D1 ; D2 = T*D2 ;
                end
        end
        
        
        function [D1,D2] = gradMat2(obj,refPoints)
        % Return the differentiation matrices so that df_dxi = Di*f(:)
        % size(Di) = [nElems nNodes] : f is defined on nodes, the gradient
        % is constant over elements
            if nargin==1 ; refPoints = obj.Points ; end
            % Build the element matrices
                Elems = [obj.Triangles NaN(size(obj.Triangles,1),1) ; obj.Quadrangles] ;
            % Coordinates
                X = refPoints(:,1) ; Y = refPoints(:,2) ;
            % Compute the mean gradient
                eee = [] ;
                nnn = [] ;
                vvvD1 = [] ;
                vvvD2 = [] ;
                for elmt = 1:size(Elems,1)
                    % Element topology
                        nodes = Elems(elmt,~isnan(Elems(elmt,:)))' ;
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
                D1 = sparse(eee(:),nnn(:),vvvD1(:),size(Elems,1),size(refPoints,1)) ;
                D2 = sparse(eee(:),nnn(:),vvvD2(:),size(Elems,1),size(refPoints,1)) ;
        end
        
        function M = tri2nod(obj)
        % Return a matrix so that f = M*F, where f is defined on Nodes and
        % F is defined on elements; size(M) = [nNodes nElems]
            M = sparse(obj.Triangles(:),reshape(repmat((1:size(obj.Triangles,1))',[1 3]),[],1),1,size(obj.Points,1),size(obj.Triangles,1)) ;
            valance = sum(M,2) ; % number of elements associated to each node
            M = M./valance(:) ;
        end
        
        function set.DataOnNodes(obj,val)
            obj.DataOnNodes = val ;
            obj.computeDataFields ;
        end
        
        function DATA = computeDataFields(obj)
        % Compute all (scalar) data fields associated to the object motion
            DATA = struct() ;
            nFrames = size(obj.MovingPoints,3) ;
            % NaNs
                DATA.NaN = obj.MovingPoints(:,1,:)*NaN ;
            % Reference Coordinates
                DATA.X1 = obj.MovingPoints(:,1,1) ;
                DATA.X2 = obj.MovingPoints(:,2,1) ;
            % Current Coordinates
                DATA.x1 = obj.MovingPoints(:,1,:) ;
                DATA.x2 = obj.MovingPoints(:,2,:) ;
            % Displacements (transformation)
                DATA.u1 = DATA.x1 - DATA.X1 ;
                DATA.u2 = DATA.x2 - DATA.X2 ;
                DATA.U = sqrt(DATA.u1.^2 + DATA.u2.^2) ;
            % Derivation matrices
                [D1,D2] = gradMat2(obj,[DATA.X1 DATA.X2]) ;
            % Transformation Gradient
                DATA.F11 = reshape(D1*squeeze(DATA.x1),[],1,nFrames) ;
                DATA.F12 = reshape(D2*squeeze(DATA.x1),[],1,nFrames) ;
                DATA.F21 = reshape(D1*squeeze(DATA.x2),[],1,nFrames) ;
                DATA.F22 = reshape(D2*squeeze(DATA.x2),[],1,nFrames) ;
            % Jacobian
                DATA.J = DATA.F11.*DATA.F22 - DATA.F12.*DATA.F21 ;
            % Cauchy Strains
                DATA.C11 = DATA.F11.*DATA.F11 + DATA.F21.*DATA.F21 ;
                DATA.C22 = DATA.F12.*DATA.F12 + DATA.F22.*DATA.F22 ;
                DATA.C12 = DATA.F11.*DATA.F12 + DATA.F21.*DATA.F22 ;
            % Euler-Lagrange Strains
                DATA.L11 = 0.5*(DATA.C11 - 1) ;
                DATA.L22 = 0.5*(DATA.C22 - 1) ;
                DATA.L12 = 0.5*(DATA.C12) ;
                DATA.L33 = -DATA.L11-DATA.L22 ;
                DATA.Lmean = 1/3*(DATA.L11+DATA.L22+DATA.L33) ;
            % Deviatoric part
                DATA.Ldev11 = DATA.L11-DATA.Lmean ;
                DATA.Ldev22 = DATA.L22-DATA.Lmean ;
                DATA.Ldev33 = DATA.L33-DATA.Lmean ;
                DATA.Ldev12 = DATA.L12 ;
            % Equivalent strain
                DATA.Leq = sqrt(2/3*(DATA.Ldev11.^2 + DATA.Ldev22.^2 + DATA.Ldev12.^2 + DATA.Ldev33.^2)) ;
            % Eigen Strains
                aaa = 1/2*(DATA.L11+DATA.L22) ;
                bbb = sqrt(1/4*(DATA.L11-DATA.L22).^2 + DATA.L12.^2) ;
                DATA.Eig1 = aaa+bbb ; % Major strain
                DATA.Eig2 = aaa-bbb ; % Minor Strain
                DATA.EigTau = bbb ; % Maximum Shear
                DATA.EigTheta = 1/2*atan(2*DATA.L12./(DATA.L11-DATA.L22)) ; % Principal Angle
            % Linearized Euler-Lagrange Strains
                DATA.E11 = reshape(D1*squeeze(DATA.u1),[],1,nFrames) ;
                DATA.E22 = reshape(D2*squeeze(DATA.u2),[],1,nFrames) ;
                DATA.E12 = 0.5*(reshape(D2*squeeze(DATA.u1),[],1,nFrames) + reshape(D1*squeeze(DATA.u2),[],1,nFrames)) ;
            % Save in the object
                obj.DataFields = DATA ;
        end
        
        
        function initSeedPreview(obj,ax)
            % INIT THE MESH
                axes(ax) ;
                triMesh = patch(obj.Points(:,1),obj.Points(:,2),NaN*obj.Points(:,2) ...
                                ,'vertices',obj.Points...
                                ,'faces',obj.Triangles...
                                ,'edgecolor','k'...
                                ,'EdgeAlpha', 0.2 ...
                                ,'linewidth',0.1...
                                ,'facecolor','interp'...
                                ,'facealpha',0.5 ...
                                ,'hittest','off' ...
                                ,'tag',obj.Name ...
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
                                letters = @(str)upper(str(uint8(str)>=65 & uint8(str)<=122)) ;
                                for f = 1:numel(fieldNames)
                                    % Add the submenu
                                        submenus(end+1) = uimenu(mData,'Label',fieldNames{f}) ;
                                    % Add a separator if needed (first letter different)
                                        if f>1
                                            if ~strcmp(letters(fieldNames{f}),letters(fieldNames{f-1}))
                                                submenus(end).Separator = 'on' ;
                                            end
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
                            submenus(end+1) = uimenu(mPlot,'Label','Points') ;
                            submenus(end+1) = uimenu(mPlot,'Label','Edges') ;
                    % Color Scale
                        mColors = uimenu(mDisplay,'Label','Colors') ;
                            mClrScale = uimenu(mColors,'Label','Scale') ;
                                submenus(end+1) = uimenu(mClrScale,'Label','Current Frame','checked','on') ;
                                submenus(end+1) = uimenu(mClrScale,'Label','All Frames') ;
                            mClrLims = uimenu(mColors,'Label','Limits') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','min-max','checked','on') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','0-max') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','symmetric') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','1*sigma') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','2*sigma') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','3*sigma') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','4*sigma') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','5*sigma') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','6*sigma') ;
                            mClrSteps = uimenu(mColors,'Label','Steps') ;
                                submenus(end+1) = uimenu(mClrSteps,'Label','Continuous','checked','on') ;
                                submenus(end+1) = uimenu(mClrSteps,'Label','11') ;
                                submenus(end+1) = uimenu(mClrSteps,'Label','7') ;
                                submenus(end+1) = uimenu(mClrSteps,'Label','4') ;
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
                    ax.UserData.axesPositionReference = [] ; % [frame posX posY]
                    ax.UserData.clrMode = 'Preset' ;
                    ax.UserData.clrScaleLabel = 'Current Frame' ;
                    ax.UserData.clrLimsLabel = 'min-max' ;
                    ax.UserData.clrStepsLabel = 'Continuous' ;
        end
        
        
        function updateSeedMenus(obj,subMenu,ax)
            % Uncheck all subMenuItems
                set(subMenu.Parent.Children,'checked','off')
            % Check the selected item
                subMenu.Checked = 'on' ;
            % Change the ax userData if needed
                switch subMenu.Parent.Label
                    case 'Data'
                        ax.UserData.dataLabel = subMenu.Label ;
                    case 'Scale'
                        ax.UserData.clrScaleLabel = subMenu.Label ;
                        ax.UserData.clrMode = 'Preset' ;
                    case 'Limits'
                        ax.UserData.clrLimsLabel = subMenu.Label ;
                        ax.UserData.clrMode = 'Preset' ;
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
                                    cen = obj.interpMat(cen,'extrap',obj.MovingPoints(:,:,frame))*obj.Points ;
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
                end
            % Update the preview
                updateSeedPreview(obj,[],ax)
        end
        
        
        function updateSeedPreview(obj,~,ax)
            % Search for the mesh
                triMesh = findobj(ax,'tag',obj.Name) ;
            % If it is not found, re-init
                if isempty(triMesh)
                    initSeedPreview(obj,ax) ;
                    triMesh = findobj(ax,'tag',obj.Name) ;
                end
            % Get the frame
                CurrentFrame = ax.UserData.currentFrame ;
            % If the frame is valid, Update
                if CurrentFrame>0 && CurrentFrame<=size(obj.MovingPoints,3)
                    % Deform the mesh
                        triMesh.Vertices = obj.MovingPoints(:,:,CurrentFrame) ;
                    % Retrieve the data
                        Data = [] ;
                        % Re-compute data fields if needed
                            if ~isfield(obj.DataFields,ax.UserData.dataLabel) ; obj.computeDataFields ; end
                        % Extract the data
                            Data = obj.DataFields.(ax.UserData.dataLabel) ;
                        % If no data, return
                            if isempty(Data) ; return ; end
                        % Re-format
                            Data = squeeze(Data) ;
                    % Apply data to mesh
                        if size(Data,1) == size(triMesh.Faces,1)
                            triMesh.FaceColor = 'flat' ;
                        else
                            triMesh.FaceColor = 'interp' ;
                        end
                        CurrentFrame = min(size(Data,2),CurrentFrame) ; % Allow a field to be constant
                        triMesh.CData = Data(:,CurrentFrame) ;
                    % COLORS
                        if strcmp(ax.UserData.clrMode,'Preset')
                            % Color scale
                                switch ax.UserData.clrScaleLabel
                                    case 'Current Frame'
                                        colorData = Data(:,CurrentFrame) ;
                                    case 'All Frames'
                                        colorData = Data ;
                                end
                                minData = min(colorData(:)) ;
                                maxData = max(colorData(:)) ;
                            % Color Limits
                                if ~any(isnan([minData maxData]))
                                    % Set CAXIS
                                        switch ax.UserData.clrLimsLabel
                                            case '0-max'
                                                if sign(minData)==sign(maxData)
                                                    if sign(minData)==-1
                                                        caxis(ax,[minData 0]) ;
                                                    else
                                                        caxis(ax,[0 maxData]) ;
                                                    end
                                                else
                                                    caxis(ax,[minData maxData]) ;
                                                end
                                            case 'min-max'
                                                caxis(ax,[minData maxData]) ;
                                            case 'symmetric'
                                                caxis(ax,max(abs([minData maxData]))*[-1 1]) ;
                                            otherwise % N*sigma
                                                % Get the factor N
                                                    N = str2double(ax.UserData.clrLimsLabel(1)) ;
                                                % Compute mean and variances
                                                    valid = ~isnan(colorData(:)) ;
                                                    avg = mean(colorData(valid)) ;
                                                    ec = std(colorData(valid)) ;
                                                % Set Color Limits
                                                    caxis(min(max(avg+N*ec*[-1 1],minData),maxData)) ;
                                        end
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
                        end
                    % AXES LIMITS
                        % Follow point if needed
                            if ~isempty(ax.UserData.axesPositionReference) && ~any(isnan(ax.UserData.axesPositionReference))
                                % Center and range of axis limits
                                    cen = round(obj.interpMat(ax.UserData.axesPositionReference(1:2))*obj.MovingPoints(:,:,ax.UserData.currentFrame)) ; 
                                    ran = round(ax.UserData.axesPositionReference(3:4)) ;
                                % Update limits
                                    [nI,nJ] = size(obj.refImgs{1}) ;
                                    ax.XLim = max(0.5,min(nJ+0.5,cen(1)+ran(1)*[-1 1]/2)) ;
                                    ax.YLim = max(0.5,min(nI+0.5,cen(2)+ran(2)*[-1 1]/2)) ;
                            end
                end
        end
        
        
    end
    
end 