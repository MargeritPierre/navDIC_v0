classdef navDICSeed_2D_Curve < navDICSeed
   
    properties
    end
    
    methods
        function obj = navDICSeed_2D_Curve(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'single') ;
                obj.Class = 'navDICSeed_2D_Curve' ;
        end
    end
    
end