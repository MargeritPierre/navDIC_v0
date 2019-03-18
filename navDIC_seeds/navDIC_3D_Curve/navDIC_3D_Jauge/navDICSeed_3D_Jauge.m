classdef navDICSeed_3D_Jauge < navDICSeed
   
    properties
        ROI = [] ;
        corrSize = [7, 7] ; % taille fenetre de correlation cam 1 et cam 2
        L0 = [] ;
    end
    
    methods
        function obj = navDICSeed_3D_Jauge(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'multiple') ;
                obj.Class = 'navDICSeed_jauge'  ;
                obj.displMethod = 'cpcorr3D' ;
                obj.strainMethod = 'deltaLCor3D' ;
           
           % Set cameras Position and Optic parameters 
               for i = 1 : length(obj.CamIDs)
                    hd = setCameraProperties(hd,i) ;
               end
           % Draw points
           nbCam = length(obj.CamIDs) ;
                for cam = 1 : nbCam
                    obj.drawToolH = drawingTool('drawROI',true ...
                                                ,'background', obj.refImgs{cam}) ;

               %  obj.drawToolH
                   if strcmpi(obj.drawToolH.Geometries(1).Class, 'impoint')
                       for p =1:2
                            obj.Points(p,:,cam) = obj.drawToolH.Geometries(p).Position ;
                       end
                   elseif strcmpi(obj.drawToolH.Geometries(1).Class, 'imline')
                       obj.Points(:,:,cam) = obj.drawToolH(cam).Geometries.Position(:,:) ;
                   end
                end
                obj.Points = camsTo3d(hd,obj.Points) ;
            % INITIALIZE
                obj.MovingPoints = ones(size(obj.Points,1),3,hd.nFrames,nbCam)*NaN ;
                obj.Displacements = ones(size(obj.Points,1),3,hd.nFrames,nbCam)*NaN ;
                obj.Strains = ones(size(obj.Points,1),1,hd.nFrames,nbCam)*NaN ;
                obj.L0 = norm(diff(obj.Points,1,1,nbCam)) ;
        end
        
        function obj = modify(obj,hd)
            for cam = 1 : length(obj.CamIDs)
                    obj.drawToolH = drawingTool('drawROI',true ...
                                                ,'background', obj.refImgs{cam}) ;

               %  obj.drawToolH
                   if strcmpi(obj.drawToolH.Geometries(1).Class, 'impoint')
                       for p =1:2
                            obj.CamPoints{cam}(p,:) = obj.drawToolH.Geometries(p).Position ;
                       end
                   elseif strcmpi(obj.drawToolH.Geometries(1).Class, 'imline')
                       obj.CamPoints{cam} = obj.drawToolH(cam).Geometries.Position(:,:) ;
                   end
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