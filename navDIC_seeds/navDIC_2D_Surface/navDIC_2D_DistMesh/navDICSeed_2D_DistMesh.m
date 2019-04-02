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
            obj.drawToolH = drawingTool(obj.drawToolH) ;
            obj.Points = obj.drawToolH.DistMesh.Points ;
            obj.Triangles = obj.drawToolH.DistMesh.Triangles ;
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
                if 0
                    clrbr = colorbar(ax) ;
                    clrbr.Units = 'pixels' ;
                    ax.Units = 'pixels' ;
                    fig = ax.Parent ;
                    margin = 0.05 ;
                    switch fig.Position(3)<fig.Position(4)
                        case true
                            clrbr.Location = 'eastoutside' ;
                            %clrbr.Position(2) = ax.Position(2) + margin*ax.Position(4) ;
                            %clrbr.Position(4) = ax.Position(4) - 2*margin*ax.Position(4) ;
                            %fig.Position(3) = fig.Position(3) + 5*clrbr.Position(3) ;
                            %clrbr.Position(1) = ax.Position(3) + clrbr.Position(3)*1 ;
                        case false
                            clrbr.Location = 'northoutside' ;
                            %fig.Position(4) = fig.Position(4) + 4*clrbr.Position(4) ;
                            %clrbr.Position(2) = ax.Position(4) + clrbr.Position(4)*1 ;
                    end
                    clrbr.AxisLocation = 'out';
                    drawnow; pause(0.05) ;
                    clrbr.Units = 'normalized' ;
                    ax.Units = 'normalized' ;
                end
            % ADD THE MENU BAR
                % Data to plot
                    mData = uimenu(ax.Parent,'text','Data') ;
                        submenus(1) = uimenu(mData,'text','|U|','checked','on') ;
                        submenus(end+1) = uimenu(mData,'text','Ux') ;
                        submenus(end+1) = uimenu(mData,'text','Uy') ;
                        submenus(end+1) = uimenu(mData,'text','Exx','separator',true) ;
                        submenus(end+1) = uimenu(mData,'text','Eyy') ;
                        submenus(end+1) = uimenu(mData,'text','Exy') ;
                % Color Scale
                    mColors = uimenu(ax.Parent,'text','Colors') ;
                        mClrScale = uimenu(mColors,'text','Scale') ;
                            submenus(end+1) = uimenu(mClrScale,'text','Current Frame','checked','on') ;
                            submenus(end+1) = uimenu(mClrScale,'text','All Frames') ;
                        mClrLims = uimenu(mColors,'text','Limits') ;
                            submenus(end+1) = uimenu(mClrLims,'text','0-max','checked','on') ;
                            submenus(end+1) = uimenu(mClrLims,'text','min-max') ;
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
                switch subMenu.Parent.Text
                    case 'Data'
                        ax.UserData.dataLabel = subMenu.Text ;
                    case 'Scale'
                        ax.UserData.clrScaleLabel = subMenu.Text ;
                    case 'Limits'
                        ax.UserData.clrLimsLabel = subMenu.Text ;
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