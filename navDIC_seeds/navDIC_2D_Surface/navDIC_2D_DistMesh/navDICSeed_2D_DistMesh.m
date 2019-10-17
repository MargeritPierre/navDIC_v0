classdef navDICSeed_2D_DistMesh < navDICSeed_2D_Surface
   
    properties
        Triangles = [] ;
        ROI = [] ;
        Shapes = [] ;
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
                obj.Displacements = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Strains = ones(size(obj.Points,1),3,hd.nFrames)*NaN ;
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
                if all(isnan(obj.MovingPoints(:))) ; 
                    MovingPoints = ones(size(Points,1),1)*obj.MovingPoints(1,:) ;
                    Displacements = MovingPoints ;
                    Strains = MovingPoints(:,[1 2 1]) ;
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
                            Displacements = ones(size(Points,1),2,hd.nFrames)*NaN ;
                            Strains = ones(size(Points,1),3,hd.nFrames)*NaN ;
                        % Project previous computations if wanted
                            if strcmp(answer,'Project')
                                % Triangulation class constructor
                                    tri = triangulation(obj.Triangles,obj.Points) ; 
                                % oldPoints associated to each newPoint
                                    % Initialization
                                        indPts = zeros(size(Points,1),3) ;
                                    % Enclosing element and barycentric positions
                                        [elmt,weights] = pointLocation(tri,Points) ; 
                                    % Outside-of-old-mesh newPoints need to be fixed (elmt=NaN)
                                        outside = isnan(elmt) ;
                                        outsidePts = Points(outside,:) ;
                                        if any(outside)
                                            % Find the closest element
                                                C = circumcenter(tri) ;
                                                [~,clstElmt] = min(sum((reshape(outsidePts,[],1,2)-reshape(C,1,[],2)).^2,3),[],2) ;
                                                elmt(outside) = clstElmt ;
                                            % Find the (extended) coordinates of the new point in each old closest element frame
                                                elmtNodes = obj.Triangles(clstElmt,:) ;
                                                p1 = obj.Points(elmtNodes(:,1),:)' ; 
                                                p2 = obj.Points(elmtNodes(:,2),:)' ; 
                                                p3 = obj.Points(elmtNodes(:,3),:)' ; 
                                                v1 = p2-p1 ; v2 = p3-p1 ; % element's frame vectors
                                                m = [] ; % coordinates in this frame
                                                for pp = 1:length(clstElmt)
                                                    m(pp,:) = [v1(:,pp) v2(:,pp)]\(outsidePts(pp,:)'-p1(:,pp)) ;
                                                end
                                            % Assign the weights
                                                indPts(outside,:) = elmtNodes ;
                                                weights(outside,:) = [1-m(:,1)-m(:,2) m(:,1) m(:,2)] ;
                                        end
                                    % Inside nodes
                                        indPts(~outside,:) = obj.Triangles(elmt(~outside),:) ;
                                % Transfer Matrix
                                    nodes = (1:size(Points,1))'*[1 1 1] ;
                                    T = sparse(nodes(:),indPts(:),weights(:),size(Points,1),size(obj.Points,1)) ;
                                % Projection of displacements
                                    MovingPoints = reshape(T*reshape(obj.MovingPoints,size(obj.Points,1),[]),size(Points,1),2,hd.nFrames) ;
                                    Displacements = reshape(T*reshape(obj.Displacements,size(obj.Points,1),[]),size(Points,1),2,hd.nFrames) ;
                                % For strain
                                    if size(obj.Strains,1)==size(obj.Points,1) % Strains defined on nodes
                                        Strains = reshape(T*reshape(obj.Strains,size(obj.Points,1),[]),size(Points,1),3,hd.nFrames) ;
                                    else % Strains defined on elements
                                        elmt(outside) = 1 ;
                                        weights(outside,:) = 1 ;
                                        Ts = sparse(1:size(Points,1),elmt,weights(:,1),size(Points,1),size(obj.Triangles,1)) ;
                                        Strains = reshape(Ts*reshape(obj.Strains,size(obj.Triangles,1),[]),size(Points,1),3,hd.nFrames) ;
                                    end
                            end
                end
                % Crush the old mesh
                    obj.drawToolH = drawToolH ;
                    obj.Points = Points ;
                    obj.Triangles = Tris ;
                    obj.MovingPoints = MovingPoints ;
                    obj.Displacements = Displacements ;
                    obj.Strains = Strains ;
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
                if 1
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
                end
            % ADD THE MENU BAR
                % Data to plot
                    mData = uimenu(ax.Parent,'Label','Data') ;
                        submenus(1) = uimenu(mData,'Label','|U|','checked','on') ;
                        submenus(end+1) = uimenu(mData,'Label','Ux') ;
                        submenus(end+1) = uimenu(mData,'Label','Uy') ;
                        submenus(end+1) = uimenu(mData,'Label','Exx','separator','on') ;
                        submenus(end+1) = uimenu(mData,'Label','Eyy') ;
                        submenus(end+1) = uimenu(mData,'Label','Exy') ;
                        submenus(end+1) = uimenu(mData,'Label','Major E','separator','on') ;
                        submenus(end+1) = uimenu(mData,'Label','Minor E') ;
                        submenus(end+1) = uimenu(mData,'Label','Max. Shear') ;
                        submenus(end+1) = uimenu(mData,'Label','Princ. Angle') ;
                % DISPLAY
                    mDisplay = uimenu(ax.Parent,'Label','Display') ;
                    % Color Scale
                        mColors = uimenu(mDisplay,'Label','Colors') ;
                            mClrScale = uimenu(mColors,'Label','Scale') ;
                                submenus(end+1) = uimenu(mClrScale,'Label','Current Frame','checked','on') ;
                                submenus(end+1) = uimenu(mClrScale,'Label','All Frames') ;
                            mClrLims = uimenu(mColors,'Label','Limits') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','0-max','checked','on') ;
                                submenus(end+1) = uimenu(mClrLims,'Label','min-max') ;
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
                % Common Properties
                    set(submenus,'callback',@(src,evt)obj.updateSeedMenus(src,ax)) ;
                % UserData in axes to choose the data to plot
                    ax.UserData.dataLabel = '|U|' ;
                    ax.UserData.clrMode = 'Preset' ;
                    ax.UserData.clrScaleLabel = 'Current Frame' ;
                    ax.UserData.clrLimsLabel = '0-max' ;
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
        
        function updateSeedPreview(obj,hd,ax)
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
                        switch ax.UserData.dataLabel
                            case '|U|'
                                Data =  sqrt(sum(obj.Displacements.^2,2)) ;
                            case 'Ux'
                                Data = obj.Displacements(:,1,:) ;
                            case 'Uy'
                                Data = obj.Displacements(:,2,:) ;
                            case 'Exx'
                                Data = obj.Strains(:,1,:) ;
                            case 'Eyy'
                                Data = obj.Strains(:,2,:) ;
                            case 'Exy'
                                Data = obj.Strains(:,3,:) ;
                            case 'Major E'
                                Data = obj.MajorStrains ;
                            case 'Minor E'
                                Data = obj.MinorStrains ;
                            case 'Max. Shear'
                                Data = obj.MaxShear ;
                            case 'Princ. Angle'
                                Data = obj.PrincipalAngle ;
                        end
                        if isempty(Data) ; return ; end
                        Data = squeeze(Data) ;
                    % Apply data to mesh
                        if size(Data,1) == size(triMesh.Faces,1) % Facet data
                            triMesh.FaceColor = 'flat' ;
                        else
                            triMesh.FaceColor = 'interp' ;
                        end
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
                end
        end
        
    end
    
end 