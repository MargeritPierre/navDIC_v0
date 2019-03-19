classdef navDICSeed_3D_Jauge < navDICSeed
   
    properties
        ROI = [] ;
        corrSize = [7, 7] ; % taille fenetre de correlation cam 1 et cam 2
        L0 = [] ;
        cam
    end
    
    methods
        function [obj] = navDICSeed_3D_Jauge(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'multiple') ;
                obj.Class = 'navDICSeed_jauge'  ;
                obj.displMethod = 'cpcorr3D' ;
                obj.strainMethod = 'deltaLCor3D' ;
           
           % Draw points
           nbCam = length(obj.CamIDs) ;
           cam = 1 ;
                while cam <= nbCam
                    obj.drawToolH = drawingTool('drawROI',true ...
                                                ,'background', obj.refImgs{cam}) ;

               %  obj.drawToolH
                   if strcmpi(obj.drawToolH.Geometries(1).Class, 'impoint')
                       for p = 1:length(obj.drawToolH.Geometries)
                            pts(p,:) = obj.drawToolH.Geometries(p).Position ;
                       end
                   elseif strcmpi(obj.drawToolH.Geometries(1).Class, 'imline')
                       pts = obj.drawToolH.Geometries.Position(:,:) ;
                   end
                   
                   if length(pts) == length(obj.Points(:,:,1)) || cam == 1 
                       obj.Points(:,:,cam) = pts ;
                       cam = cam + 1 ;
                   end
                   
                end
            % INITIALIZE
                obj.MovingPoints = ones(size(obj.Points,1),size(obj.Points,2),hd.nFrames,nbCam)*NaN ;
                obj.Displacements = ones(size(obj.Points,1),size(obj.Points,2),hd.nFrames,nbCam)*NaN ;
                obj.Strains = ones(1,1,hd.nFrames,nbCam)*NaN ;
                for i = 1:nbCam
                    obj.L0(i) = norm(diff(obj.Points(:,:,i),1,1)) ;
                end
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
                    idCam = hd.Previews{1}.cam ;
                    line.XData = obj.MovingPoints(:,1,hd.CurrentFrame,idCam) ;
                    line.YData = obj.MovingPoints(:,2,hd.CurrentFrame,idCam) ;
                    label.Position = mean(obj.MovingPoints(:,:,hd.CurrentFrame,idCam),1) ;
                    label.String = [num2str(obj.Strains(1,1,hd.CurrentFrame,idCam)*100,3), '%'] ;
                %end
            end
        end
        
    end
    
end