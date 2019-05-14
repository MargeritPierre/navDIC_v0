classdef navDICSeed_2D_Line < navDICSeed_2D_Curve
   
    properties
    end
    
    methods
        function obj = navDICSeed_2D_Line(hd)
            % Initialize
                obj = obj@navDICSeed_2D_Curve(hd) ;
                obj.Class = 'navDICSeed_2D_Line' ;
                obj.strainMethod = 'lineStrain' ;
           % Draw points       
                uiContextMenu = {'N',num2cell(num2str(round((2).^(1:12)')),2),{'16'} ...
                                            ,'CORRSIZEX',num2cell(num2str((5:40)'),2),{'10'} ...
                                            ,'CORRSIZEY',num2cell(num2str((5:40)'),2),{'10'} ...
                                            } ;
                
                obj.drawToolH = drawingTool('drawROI',true ...
                                            ,'background', obj.refImgs{1} ...
                                            ,'uicontextmenu',uiContextMenu ...
                                            ,'updateCallback',@(H)navDIC_processShapesForLines(H)...
                                            ) ;
                obj.Points = obj.drawToolH.Points ;
           % Process drawing tool output
               obj.Points = obj.drawToolH.Points ;
            % INITIALIZE
                obj.MovingPoints = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Displacements = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Strains = ones(size(obj.Points,1)-1,1,hd.nFrames)*NaN ;
        end
        
        function obj = modify(obj,hd)
            obj.drawToolH = drawingTool(obj.drawToolH) ;
            obj.Points = obj.drawToolH.Points ;
            % INITIALIZE
                obj.MovingPoints = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Displacements = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Strains = ones(size(obj.Points,1)-1,1,hd.nFrames)*NaN ;
        end
        
        function updateSeedPreview(obj,hd,ax)
            pts = findobj(ax,'tag',obj.Name) ;
            if isempty(pts)
                axes(ax) ;
                pts = plot(obj.Points(:,1),obj.Points(:,2) ...
                                ,'.-b'...
                                ,'markersize', 20 ...
                                ,'linewidth',1 ...
                                ,'hittest','off' ...
                                ,'tag',obj.Name ...
                                ) ;
            end
            if hd.CurrentFrame>0
                n = length(obj.MovingPoints(:,1,hd.CurrentFrame)) ;
                pts.XData = obj.MovingPoints(:,1,hd.CurrentFrame) ;
                pts.YData = obj.MovingPoints(:,2,hd.CurrentFrame) ;
                %cd = single(obj.Strains(:,1,hd.CurrentFrame)) ;
                %drawnow
                %set(pts.Edge, 'ColorBinding','interpolate', 'ColorData',cd)
            end
        end
        
    end
    
end