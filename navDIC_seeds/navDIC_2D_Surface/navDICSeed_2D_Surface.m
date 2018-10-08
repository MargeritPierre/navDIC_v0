classdef navDICSeed_2D_Surface < navDICSeed
   
    properties
        corrSize = 40 ;
    end
    
    methods
        function obj = navDICSeed_2D_Surface(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'single') ;
                obj.Class = 'navDICSeed_2D_Surface' ;
        end
        function obj = modify(obj,hd)
        end
        function updateSeedPreview(obj,ax)
        end
    end
    
end