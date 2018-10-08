classdef navDICSeed_1D_Jauge < navDICSeed
   
    properties
    end
    
    methods
        function obj = navDICSeed_1D_Jauge(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'single') ;
                obj.Class = 'navDICSeed_jauge' ;

                
           % Draw points       
                obj.drawToolH = drawingTool('drawROI',true ...
                                            ,'background', obj.refImgs{1}) ;
           %  obj.drawToolH
                if strcmp(obj.drawToolH.Geometry.Class, 'imline')
                    obj.Points(:,:) = obj.drawToolH.Geometries.Positions(:,:) ;
                elseif strcmp(obj.drawToolH.Geometry.Class, 'impoint')
                    for p =1:2
                        obj.Points(p,:) = obj.drawToolH.Geometries(p).Position ;
                    end
                end
           obj.strainMethod = 'deltaL' ;
           obj.displacementMethod = 'cpcorr' ;
        end
        function obj = modify(obj,hd)
            obj.drawToolH = drawingTool(obj.drawToolH) ;
            if strcmp(obj.drawToolH.Geometry.Class, 'imline')
                obj.Points(:,:) = obj.drawToolH.Geometries.Positions(:,:) ;
            elseif strcmp(obj.drawToolH.Geometry.Class, 'impoint')
                for p =1:2
                    obj.Points = obj.drawToolH.Geometries(p).Positions ;
                end
            end
            
        end
    end
    
end