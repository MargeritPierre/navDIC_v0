classdef navDICPreview
    
    properties
        fig = [] ; % the main figure
        handles = [] ; % handles linked to the preview
        isValid = true ; 
    end
    
    methods(Static)
        
        % CONSTRUCTOR
            function prev = navDICPreview(hd)
                % Parameters
                    figRelSize = 0.5 ;
                % Figure Creation
                    prev.fig = figure('tag',hd.navDICTag) ;
                    prev.fig.ToolBar = 'none' ;
                    prev.fig.MenuBar = 'none' ;
                    prev.fig.NumberTitle = 'off' ;
                    prev.fig.Name = 'navDIC Default Preview' ;
                % Positionning
                    prev.fig.Units = 'normalized' ;
                    prev.fig.Position = [.5-figRelSize/2 .5-figRelSize/2 figRelSize figRelSize] ;
            end
        
        % UPDATE
            function prev = updatePreview(prev,hd)
                % Test the preview validity
                    if ~isvalid(prev.fig)
                        prev.isValid = false ;
                        return ;
                    end
            end
        
        % DESTRUCTOR
            function hd = closePreview(prev,hd)

            end
    end
    
end