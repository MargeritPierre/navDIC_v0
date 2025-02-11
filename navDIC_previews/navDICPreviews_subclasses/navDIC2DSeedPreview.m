classdef navDIC2DSeedPreview < navDICCameraPreview
    
    properties
        Seed
    end
    
    methods
        
        % CONSTRUCTOR        
            function prev = navDIC2DSeedPreview(hd)
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
                    prev = prev@navDICCameraPreview(hd,camID,false) ;
                    if ~prev.isValid ; return ; end
                    prev.Seed = seed ;
                    prev.fig.Name = ['navDIC Seed Preview: '+string(prev.Seed.Name)] ;
                % Update the preview
                    prev = updatePreview(prev,hd) ;
            end
            
        % UPDATE
            function [prev,hd] = updatePreview(prev,hd)
                % Superclass updating function
                    prev = updatePreview@navDICCameraPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
                % Update seed preview
                    prev.AxesImg.UserData.currentFrame = hd.CurrentFrame ;
                    prev.Seed.updateSeedPreview(hd,prev.AxesImg) ;
            end
    end
    
            
    % OTHER FUNCTIONS
    
        methods(Static)
                
        end
    
end