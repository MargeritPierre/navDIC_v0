classdef navDICSeed_1D_Jauge < navDICSeed
   
    properties
        ROI = [] ;
        corrSize = 40 ;
        L0 = [] ;
    end
    
    methods
        
        function obj = navDICSeed_1D_Jauge(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'single') ;
                obj.Class = 'navDICSeed_jauge' ;
                obj.strainMethod = 'deltaL' ;
                
           % Draw points       
                obj.drawToolH = drawingTool('drawROI',true ...
                                            ,'background', obj.refImgs{1}) ;
           %  obj.drawToolH
               if strcmpi(obj.drawToolH.Geometries.Class, 'impoint')
                   for p =1:2
                        obj.Points(p,:) = obj.drawToolH.Geometries(p).Position ;
                   end
               elseif strcmpi(obj.drawToolH.Geometries.Class, 'imline')
                   obj.Points(:,:) = obj.drawToolH.Geometries.Position(:,:) ;
               end
            % INITIALIZE
                obj.MovingPoints = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Displacements = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Strains = ones(size(obj.Points,1),1,hd.nFrames)*NaN ;
                obj.L0 = norm(diff(obj.Points,1,1)) ;
        end
        
        function obj = modify(obj,hd)
            obj.drawToolH = drawingTool(obj.drawToolH) ;
            %  obj.drawToolH
           if strcmpi(obj.drawToolH.Geometries.Class, 'impoint')
               for p =1:2
                    obj.Points(p,:) = obj.drawToolH.Geometries(p).Position ;
               end
           elseif strcmpi(obj.drawToolH.Geometries.Class, 'imline')
               obj.Points(:,:) = obj.drawToolH.Geometries.Position(:,:) ;
           end
        end
        
        function obj = updateSeedPreview(obj,hd,ax)
            handles = findobj(ax,'tag',obj.Name) ;
            line = findobj(handles,'type','line') ;
            label = findobj(handles,'type','text') ;
            if isempty(line)
                line = plot(ax ...
                                ,NaN,NaN ...
                                ,'linewidth',2 ...
                                ,'marker','o'...
                                ,'color','b'...
                                ,'hittest','off' ...
                                ,'tag',obj.Name ...
                                ) ;
                axes(ax) ;
                label = text(...ax, ...
                                NaN,NaN,'' ...
                                ,'color',[0 0 1]+[1 1 0]*.5 ...
                                ,'fontsize',13 ...
                                ,'backgroundcolor','w' ...
                                ,'interpreter','none' ...
                                ,'horizontalalignment','left'...
                                ,'hittest','off' ...
                                ,'tag',obj.Name ...
                                ) ;
            end
            if hd.CurrentFrame>0
                %try
                    line.XData = obj.MovingPoints(:,1,hd.CurrentFrame) ;
                    line.YData = obj.MovingPoints(:,2,hd.CurrentFrame) ;
                    label.Position = mean(obj.MovingPoints(:,:,hd.CurrentFrame),1) ;
                    label.String = [num2str(obj.Strains(1,1,hd.CurrentFrame)*100,3), '%'] ;
                %end
            end
        end
        
    end
    
end