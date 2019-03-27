 classdef navDICSeed_3D_Jauge < navDICSeed
   
    properties
        ROI = [] ;
        corrSize = [6, 6] ; % taille fenetre de correlation cam 1 et cam 2
        L0 = [] ;
        cam = [] ;
        t = [] ;
    end
    
    methods
        function [obj] = navDICSeed_3D_Jauge(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'multiple') ;
                obj.Class = 'navDICSeed_jauge'  ;
                obj.displMethod = 'cpcorr3D' ;
                obj.strainMethod = 'deltaLCor3DRT' ;
           
                
           % Draw points
           nbCam = length(obj.CamIDs) ;
           cam = 1 ;
                while cam <= nbCam
%                     refpts = [] ;
%                     while isempty(refpts)
%                         obj.drawToolH = drawingTool('drawROI',true ...
%                                                 ,'background', obj.refImgs{cam},'title',['DrawingTool : Placer le points de reference de la camera :', num2str(cam)]) ;
%                         refpts = obj.drawToolH.Geometries(1).Position ;
%                     end
                    obj.drawToolH = drawingTool('drawROI',true ...
                                                ,'background', obj.refImgs{cam},'title',['DrawingTool : Placer les jauges de la cam?ra :', num2str(cam)],'corrSize','corrSize') ;
                    pts = [] ;
                   
               %  obj.drawToolH
                   if strcmpi(obj.drawToolH.Geometries(1).Class, 'impoint')
                       for p = 1:length(obj.drawToolH.Geometries)
                            pts(p,:) = obj.drawToolH.Geometries(p).Position ;
                       end
                   elseif strcmpi(obj.drawToolH.Geometries(1).Class, 'imline')
                       for l = 1:length(obj.drawToolH.Geometries)
                            pts(end+1:end+2,:) = obj.drawToolH.Geometries(l).Position(:,:) ;
                       end
                   end
                   
                   if cam == 1 || length(pts) == length(obj.Points(:,:,1))
                       obj.Points(:,:,cam) = pts ;
                       cam = cam + 1 ;
                   end
                   
                end
            % INITIALIZE
                obj.MovingPoints = ones(size(obj.Points,1),size(obj.Points,2),hd.nFrames,nbCam)*NaN ;
                obj.Displacements = ones(size(obj.Points,1),size(obj.Points,2),hd.nFrames,nbCam)*NaN ;
                obj.Strains = ones(size(obj.Points,1)/2,2,hd.nFrames,nbCam)*NaN ;
                
                refPts = camsTo3d(hd,obj.Points) ;

                for i = 1:nbCam
                    obj.L0(:,i) = sqrt(sum((refPts(2:2:end,:,i) - refPts(1:2:end,:,i)).^2,2)) ;
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
                for i = 1:length(line)
                    camName = strsplit(line(i).Parent.Parent.Name, 'Camera ') ;
                    camName = camName{end} ;
                    idCam = find(ismember({hd.Cameras.Name},camName)) ;
                    line(i).XData = obj.MovingPoints((i-1)*2+1:i*2,1,hd.CurrentFrame,idCam) ;
                    line(i).YData = obj.MovingPoints((i-1)*2+1:i*2,2,hd.CurrentFrame,idCam) ;
                    label(i).Position = mean(obj.MovingPoints((i-1)*2+1:i*2,:,hd.CurrentFrame,idCam),1) ;
                    label(i).String = [num2str(obj.Strains(i,1,hd.CurrentFrame,idCam)*100,3), '%'] ;
                end
                
            end
        end
        
    end
    
end