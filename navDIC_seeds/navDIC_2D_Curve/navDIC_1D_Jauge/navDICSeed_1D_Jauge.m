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
           for p =1:2
                obj.Points(p,:) = obj.drawToolH.Geometries(p).Position ;
           end
           obj.strainMethod = 'deltaL_L' ;
        end
        function obj = modify(obj,hd)
                obj.drawToolH = drawingTool(obj.drawToolH) ;
                obj.Points = obj.drawToolH.Positions ;
        end
        
        
    end
    
end