classdef navDIC2DSeedPreview < navDICCameraPreview
    
    properties
        SeedName = {} ;
        cam = [] ;
    end
    
    methods
        
        % CONSTRUCTOR        
            function prev = navDIC2DSeedPreview(hd)
                % Choose a seed to preview
                    [seedID,valid] = selectSeeds(hd,'single') ;
                % get the cam to preview
                    if ~valid 
                        camID = -1 ;
                    else
                        seed = hd.Seeds(seedID) ;
                        camID = seed.CamIDs ;
                    end
                    if length(camID) > 1
                        listCam = {hd.Cameras(:).Name} ;
                        ID = listdlg('PromptString','Select the Camera :',...
                                'SelectionMode','Single',...
                                'initialValue',1,...
                                'ListString',listCam) ;
                    else
                        ID = camID ; 
                    end
                % Superclass constructor call
                    prev = prev@navDICCameraPreview(hd,ID) ;
                    prev.fig.Name = ['navDIC Seed Preview: ',seed.Name, ', Camera ', hd.Cameras(ID).Name] ;
                    prev.SeedName = seed.Name ;
                    prev.cam = ID ;
                    [prev,hd] = setIdPreview(prev,hd) ;

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
            
            function [prev,hd] = setIdPreview(prev,hd)
                % Superclass updating function
                    seedID = ismember({hd.Seeds.Name},prev.SeedName) ;
                    hd.Seeds(seedID).cam = prev.cam ;
            end
            
    end
    
            
    % OTHER FUNCTIONS
    
        methods(Static)
                
        end
    
end