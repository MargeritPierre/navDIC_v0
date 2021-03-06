classdef navDIC2DSeedGraph < navDICCameraPreview
    
    properties
        SeedName = {} ;
    end
    
    methods
        
        % CONSTRUCTOR        
            function prev = navDIC2DSeedGraph(hd)
                % Choose a seed to preview
                    [seedID,valid] = selectSeeds(hd,'single') ;
                % Get the cam to preview
                    if ~valid 
                        camID = -1 ;
                    else
                        seed = hd.Seeds(seedID) ;
                        camID = seed.CamIDs ;
                    end
                % Superclass constructor call
                    prev = prev@navDICCameraPreview(hd,camID) ;
                    prev.fig.Name = ['navDIC Seed Preview: ',seed.Name] ;
                    prev.SeedName = seed.Name ;
                % Update the preview
                    prev = updatePreview(prev,hd) ;
            end
            
        % UPDATE
            function prev = updatePreview(prev,hd)
                % Superclass updating function
                    prev = updatePreview@navDICCameraPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
                % Update seed preview
                    seedID = ismember({hd.Seeds.Name},prev.SeedName) ;
                    hd.Seeds(seedID).updateSeedPreview(hd,prev.handles.AxesImg) ;
            end
    end
    
            
    % OTHER FUNCTIONS
    
        methods(Static)
                
        end
    
end