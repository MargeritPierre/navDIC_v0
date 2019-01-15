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
        
        function updateSeedPreview(obj,hd,ax)
            triMesh = findobj(ax,'tag',obj.Name) ;
            if isempty(triMesh)
                axes(ax) ;
                triMesh = patch(obj.Points(:,1),obj.Points(:,2),NaN*obj.Points(:,2) ...
                                ,'vertices',obj.Points...
                                ,'faces',obj.Triangles...
                                ,'edgecolor','k'...
                                ,'EdgeAlpha', 0.2 ...
                                ,'linewidth',0.1...
                                ,'facecolor','interp'...
                                ,'facealpha',0.7 ...
                                ,'hittest','off' ...
                                ,'tag',obj.Name ...
                                ) ;
                % ADD A COLORBAR
                    if 1
                        clrbr = colorbar(ax,'location','west') ;
                        clrbr.Units = 'pixels' ;
                        ax.Units = 'pixels' ;
                        fig = ax.Parent ;
                        fig.Position(3) = fig.Position(3) + 5*clrbr.Position(3) ;
                        clrbr.Position(1) = ax.Position(3) + clrbr.Position(3)*1 ;
                        clrbr.AxisLocation = 'out';
                        drawnow;
                        clrbr.Units = 'normalized' ;
                        ax.Units = 'normalized' ;
                    end
            end
            if hd.CurrentFrame>0
                triMesh.Vertices = obj.MovingPoints(:,:,hd.CurrentFrame) ;
                Data =  ... obj.Strains(:,2,:) ... Eyy
                        ... obj.Strains(:,1,:) ... Exx
                         obj.Strains(:,3,:) ... Exy
                        ... obj.Displacements(:,1,:) ... Ux
                        ... obj.Displacements(:,2,:) ... Uy
                        ... sqrt(sum(obj.Displacements(:,:,:).^2,2)) ... Displacement Magnitude
                        ;
                %Data = log10(abs(Data)+1) ;
                Data = squeeze(Data) ;
                triMesh.CData = Data(:,hd.CurrentFrame) ;
                % Color scale
                    minData = min(Data(:)) ;
                    maxData = max(Data(:)) ;
                    if sign(minData)==sign(maxData) ;
                        if sign(minData)==-1
                            caxis(ax,[min(Data(:)) 0]) ;
                        else
                            caxis(ax,[0 max(Data(:))]) ;
                        end
                    else
                        caxis(ax,[min(Data(:)) max(Data(:))]) ;
                    end
            end
        end
        
    end
    
end 