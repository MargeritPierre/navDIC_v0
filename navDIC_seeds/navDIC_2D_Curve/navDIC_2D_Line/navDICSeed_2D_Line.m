classdef navDICSeed_2D_Line < navDICSeed_2D_Curve
   
    properties
    end
    
    methods
        function obj = navDICSeed_2D_Line(hd)
            % Initialize
                obj = obj@navDICSeed_2D_Curve(hd) ;
                obj.Class = 'navDICSeed_2D_Line' ;
        end
    end
    
end