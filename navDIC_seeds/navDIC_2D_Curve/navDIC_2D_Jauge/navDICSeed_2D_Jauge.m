classdef navDICSeed_2D_Jauge < navDICSeed
   
    properties
        ROI = [] ;
        corrSize = 8 ;
        L0 = [] ;
        A = [] ;
    end
    
    
    methods
        
        function obj = navDICSeed_2D_Jauge(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'single') ;
                obj.Class = 'navDICSeed_jauge' ;
                obj.strainMethod = 'deltaL' ;
                
           % Draw points       
                obj.drawToolH = drawingTool('drawROI',true ...
                                            ,'background', obj.refImgs{1}) ;
           %  obj.drawToolH
               
               if strcmpi(obj.drawToolH.Geometries(1).Class, 'impoint')
                   for p = 1:length(obj.drawToolH.Geometries)
                        obj.Points(p,:) = obj.drawToolH.Geometries(p).Position ;
                   end
               elseif strcmpi(obj.drawToolH.Geometries(1).Class, 'imline')
                   for l = 1:length(obj.drawToolH.Geometries)
                        obj.Points(end+1:end+2,:) = obj.drawToolH.Geometries(l).Position(:,:) ;
                   end
               end
            % INITIALIZE
                obj.MovingPoints = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Displacements = ones(size(obj.Points,1),2,hd.nFrames)*NaN ;
                obj.Strains = ones(size(obj.Points,1)/2,1,hd.nFrames)*NaN ;
                obj.L0 = norm( obj.Points(2:2:end,:) - obj.Points(1:2:end,:) ) ;
        end
        
        function obj = modify(obj,hd)
            obj.drawToolH = drawingTool(obj.drawToolH) ;
            %  obj.drawToolH
           if strcmpi(obj.drawToolH.Geometries(1).Class, 'impoint')
               for p = 1:length(obj.drawToolH.Geometries)
                   obj.Points(p,:) = obj.drawToolH.Geometries(p).Position ;
               end
           elseif strcmpi(obj.drawToolH.Geometries(1).Class, 'imline')
               for l = 1:length(obj.drawToolH.Geometries)
                   obj.Points(end+1:end+2,:) = obj.drawToolH.Geometries(l).Position(:,:) ;
               end
           end
        end
        
        function obj = updateSeedPreview(obj,hd,ax)
            handles = findobj(ax,'tag',obj.Name) ;
            line = findobj(handles,'type','line') ;
            label = findobj(handles,'type','text') ;
            if isempty(line)
                axes(ax) ;
                for l = 1:size(obj.Strains,1)
                    line(l) = plot(ax ...
                                    ,NaN,NaN ...
                                    ,'linewidth',2 ...
                                    ,'marker','o'...
                                    ,'color','b'...
                                    ,'hittest','off' ...
                                    ,'tag',obj.Name ...
                                    ) ;
                    
                    label(l) = text(...ax, ...
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
            end
            if hd.CurrentFrame>0
                %try
                for i = 1:length(line)
                    line(i).XData = obj.MovingPoints((i-1)*2+1:i*2,1,hd.CurrentFrame) ;
                    line(i).YData = obj.MovingPoints((i-1)*2+1:i*2,2,hd.CurrentFrame) ;
                    label(i).Position = mean(obj.MovingPoints((i-1)*2+1:i*2,:,hd.CurrentFrame),1) ;
                    label(i).String = [num2str(obj.Strains(i,1,hd.CurrentFrame)*100,3), '%'] ;
                end
                %end
            end
        end
        
    end
    
end