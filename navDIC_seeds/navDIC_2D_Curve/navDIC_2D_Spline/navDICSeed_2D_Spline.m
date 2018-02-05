classdef navDICSeed_2D_Spline < navDICSeed_2D_Curve
   
    properties
    end
    
    methods
        function obj = navDICSeed_2D_Spline(hd)
            % Initialize
                obj = navDICSeed_2D_Curve(hd) ;
                obj.Class = 'navDICSeed_2D_Spline' ;
        end
    end
    
end