classdef navDICPlotPreview < navDICPreview
    
    properties
    end
    
    methods
        
        % CONSTRUCTOR        
            function prev = navDICPlotPreview(hd,varargin)
                % Superclass constructor call
                    prev = prev@navDICPreview(hd,varargin{:}) ;
                    prev.fig.Name = 'navDIC Plot Preview' ;
                % Set the axes
                    prev.Axes = axes('outerposition',[0 0 1 1]) ;
            end
            
        % UPDATE
            function prev = updatePreview(prev,hd)
                % Superclass updating function
                    prev = updatePreview@navDICPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
            end
        
        % DESTRUCTOR
            function closePreview(obj)
            end
            
    end
    
            
    % OTHER FUNCTIONS
    
        methods(Static)                
        end
    
end