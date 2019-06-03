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
                uiContextMenu = {'h0',num2cell(num2str(round((1.4).^(1:15)')),2),{'29'};...
                    'h',num2cell(num2str(round((1.4).^(1:15)')),2),{'29'};...
                    'l',num2cell(num2str(round((1.4).^(8:23)')),2),{'79'};...
                    } ;
                obj.drawToolH = drawingTool('drawROI',true ...
                                            ,'background', obj.refImgs{1} ...
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
            % ASK THE USER WHAT TO DO
                if all(isnan(obj.MovingPoints(:))) ; return ; end
                answer = questdlg(...
                            {'THE MESH HAS BEEN MODIFIED.',...
                            'WHAT DO YOU WANT TO DO WITH THE EXISTING DATA ?'}...
                            ,'/!\ WARNING'...
                            ,'Project','Clear','Cancel','Project'...
                            ) ; 
                if strcmp(answer,'Cancel') ; return ; end
            % MODIFY THE MESH DATA
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
                                outsideNodes = isnan(elmt) ;
                                % Assign the nearest old vertex
                                    indPts(outsideNodes,:) = repmat(nearestNeighbor(tri,Points(outsideNodes,:)),[1 3]) ;
                                    weights(outsideNodes,:) = 1/3 ;
                            % Inside nodes
                                indPts(~outsideNodes,:) = obj.Triangles(elmt(~outsideNodes),:) ;
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
                                elmt(outsideNodes) = 1 ;
                                weights(outsideNodes,:) = NaN ;
                                Ts = sparse(1:size(Points,1),elmt,weights(:,1),size(Points,1),size(obj.Triangles,1)) ;
                                Strains = reshape(Ts*reshape(obj.Strains,size(obj.Triangles,1),[]),size(Points,1),3,hd.nFrames) ;
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
                    %clrbr.Units = 'pixels' ;
                    %ax.Units = 'pixels' ;
                    fig = ax.Parent ;
                    %margin = 0.05 ;
                    switch fig.Position(3)<fig.Position(4)
                        case true
                            clrbr.Location = 'east' ;
                            %clrbr.Position(2) = ax.Position(2) + margin*ax.Position(4) ;
                            %clrbr.Position(4) = ax.Position(4) - 2*margin*ax.Position(4) ;
                            %fig.Position(3) = fig.Position(3) + 5*clrbr.Position(3) ;
                            %clrbr.Position(1) = ax.Position(3) + clrbr.Position(3)*1 ;
                        case false
                            clrbr.Location = 'north' ;
                            %fig.Position(4) = fig.Position(4) + 4*clrbr.Position(4) ;
                            %clrbr.Position(2) = ax.Position(4) + clrbr.Position(4)*1 ;
                    end
                    %clrbr.AxisLocation = 'out';
                    %drawnow; pause(0.05) ;
                    %clrbr.Units = 'normalized' ;
                    %ax.Units = 'normalized' ;
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
                % Color Scale
                    mColors = uimenu(ax.Parent,'Label','Colors') ;
                        mClrScale = uimenu(mColors,'Label','Scale') ;
                            submenus(end+1) = uimenu(mClrScale,'Label','Current Frame','checked','on') ;
                            submenus(end+1) = uimenu(mClrScale,'Label','All Frames') ;
                        mClrLims = uimenu(mColors,'Label','Limits') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','0-max','checked','on') ;
                            submenus(end+1) = uimenu(mClrLims,'Label','min-max') ;
                % Common Properties
                    set(submenus,'callback',@(src,evt)obj.updateSeedMenus(src,ax)) ;
                % UserData in axes to choose the data to plot
                    ax.UserData.dataLabel = '|U|' ;
                    ax.UserData.clrScaleLabel = 'Current Frame' ;
                    ax.UserData.clrLimsLabel = '0-max' ;
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
                    case 'Limits'
                        ax.UserData.clrLimsLabel = subMenu.Label ;
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
                        switch ax.UserData.dataLabel
                            case '|U|'
                                Data =  sqrt(sum(obj.Displacements(:,:,:).^2,2)) ;
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
                        end
                        Data = squeeze(Data) ;
                    % Apply data to mesh
                        if size(Data,1) == size(triMesh.Faces,1) % Facet data
                            triMesh.FaceColor = 'flat' ;
                        else
                            triMesh.FaceColor = 'interp' ;
                        end
                        triMesh.CData = Data(:,CurrentFrame) ;
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
                            if sign(minData)==sign(maxData)
                                if sign(minData)==-1
                                    caxis(ax,[minData 0]) ;
                                else
                                    caxis(ax,[0 maxData]) ;
                                end
                            else
                                caxis(ax,[minData maxData]) ;
                            end
                        end
                end
        end
        
    end
    
end 