classdef navDICSeed
    
    properties
        % Basic Infos
            Name = 'defaultSeedName' ;
            Class = 'navDICSeed' ;
            isValid = true ;
        % Image sources
            CamIDs = [] ;
            refImgs = {} ;
        % Seed elements
            Points = [] ;
            MovingPoints = [] ;
            Displacements = [] ;
            Strains = [] ;
        % Displ. Computation Method
            displMode = 'rel' ;
            displMethod = 'cpcorr' ;
        % Strain Computation Method
            strainMethod = 'planefit' ;
        % drawingTool retrurns
            drawToolH = [] ;
    end
    
    methods
        function obj = navDICSeed(hd,nCams)
            if ~isempty(hd.Cameras) % Is there defined cameras ?
                % Select related Cameras
                    [IDs,validIDs] = selectCameras(hd,nCams) ;
                    obj.CamIDs = IDs ;
                    obj.isValid = obj.isValid && validIDs ;
                    if ~validIDs ; return ; end
                % Get reference Images
                    for id = IDs
                        obj.refImgs{end+1} = im2single(getsnapshot(hd.Cameras(id).VidObj)) ;
                    end
            elseif hd.debug % Debug mode
                % Simulate a camera
                    vidRes = [640 480] ;
                    nCams = 1 ;
                    for id = 1:nCams
                        obj.refImgs{end+1} = ones(vidRes)'*0.5 ;
                    end
            else % The seed is Invalid
                obj.isValid = false ;
            end
        end
        function obj = modify(obj,hd)
        end
        function obj = computeDisplacements(obj,hd)
            script = ['navDIC_',obj.displMethod] ;
            eval(['[obj,hd] = ',script,'(obj,hd) ;']) ;
        end
        function obj = computeStrains(obj,hd)
            script = ['navDIC_',obj.strainMethod] ;
            eval(['[obj,hd] = ',script,'(obj,hd) ;']) ;
        end
        function updateSeedPreview(obj,ax)
        end
    end
    
end