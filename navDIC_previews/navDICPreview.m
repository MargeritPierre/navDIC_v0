classdef navDICPreview
    
    properties
        fig = [] ; % the main figure
        isValid = true ; 
    end
    
    methods
        
        % CONSTRUCTOR
            function prev = navDICPreview(hd,varargin)
                % Parameters
                    figRelSize = 0.7 ; % relative size to screen
                % Figure Creation
                    prev.fig = figure('tag',hd.navDICTag) ;
                    prev.fig.MenuBar = 'none' ; % 'figure' ;
                    prev.fig.ToolBar = 'figure' ; % 'none' ;
                    prev.fig.NumberTitle = 'off' ;
                    prev.fig.Name = 'navDIC Default Preview' ;
                % Positionning
                    %prev.fig.Units = 'normalized' ;
                    L = prev.fig.Position(3:4) ;
                    prev.fig.Position = prev.fig.Position + (1-figRelSize)/2*[L 0 0] - (1-figRelSize)*[0 0 L] ;
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