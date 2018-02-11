classdef navDICPreview
    
    properties
        fig = [] ; % the main figure
        isValid = true ; 
    end
    
    methods
        
        % CONSTRUCTOR
            function prev = navDICPreview(hd,varargin)
                % Parameters
                    figRelSize = 0.5 ; % relative size to screen
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
            function prev = updatePreview(prev,hd,varargin)
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