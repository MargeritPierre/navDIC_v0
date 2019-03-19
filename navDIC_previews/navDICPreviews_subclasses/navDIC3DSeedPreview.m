classdef navDIC3DSeedPreview < navDICCameraPreview
    
    properties
        SeedName = {} ;
    end
    
    methods
        
        % CONSTRUCTOR        
            function prev = navDIC3DSeedPreview(hd)
                % Choose a seed to preview
                    [seedID,valid] = selectSeeds(hd,'single') ;
                % get the cam to preview
                    if ~valid 
                        camID = -1 ;
                    else
                        seed = hd.Seeds(seedID) ;
                        camID = seed.CamIDs ;
                    end
                % Superclass constructor call
                    prev = prev@navDICCameraPreview(hd,camID) ;
                    prev.fig.Name = ['navDIC Seed Preview: ',seed.Name, ' camera : ', num2str(camID(i)) ] ;
                    prev.SeedName = seed.Name ;
            end
            
        % UPDATE
            function [prev,hd] = updatePreview(prev,hd)
                % Superclass updating function
                    prev = updatePreview@navDICCameraPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
                % Update seed preview
                    seedID = ismember({hd.Seeds.Name},prev.SeedName) ;
                    hd.Seeds(seedID).updateSeedPreview(hd,prev.AxesImg) ;
            end
    end
    
            
    % OTHER FUNCTIONS
    
        methods(Static)
                
        end
    
end