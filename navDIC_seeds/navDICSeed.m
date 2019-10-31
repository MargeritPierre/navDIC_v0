classdef navDICSeed < handle & matlab.mixin.Heterogeneous
    
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
            MovingPoints_bkp = [] ;
            Displacements = [] ;
            Strains = [] ;
            MajorStrains = [] ;
            MinorStrains = [] ;
            MaxShear = [] ;
            PrincipalAngle = [] ;
        % Displ. Computation Method
            compDisp = true ;
            displMode = 'abs' ;
            useExistingDisp = true ;
            RefFrame = 1 ;
            displMethod = ... 'cpcorr' ...
                           'fftdisp' ...
                          ;
        % Strain Computation Method
            compStrains = true ;
            strainMethod = 'planefit' ;
        % drawingTool retrurns
            drawToolH = [] ;
    end
    
    methods
        
        function obj = navDICSeed(hd,nCams)
            % CAMERA-DEPENDENT DATA
                if ~isempty(hd.Cameras) % Is there defined cameras ?
                    % Select related Cameras
                        [IDs,validIDs] = selectCameras(hd,nCams) ;
                        obj.CamIDs = IDs ;
                        obj.isValid = obj.isValid && validIDs ;
                        if ~validIDs ; return ; end
                    % Get reference Images
                        for id = IDs
                            if ~isempty(hd.Images)
                                obj.refImgs{end+1} = hd.Images{id}(:,:,:,hd.CurrentFrame) ;
                            else
                                hd = startAllCameras(hd) ;
                                obj.refImgs{end+1} = im2single(getsnapshot(hd.Cameras(id).VidObj)) ;
                            end
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
            % SEED PARAMETERS
                obj.RefFrame = max(1,hd.CurrentFrame) ;
                % If there was token frames before this seed creation
                    if hd.CurrentFrame>1
                        obj.MovingPoints = repmat(obj.Points,[1 1 obj.RefFrame]) ;
                        obj.Displacements = zeros([size(obj.Points) obj.RefFrame]) ;
                        obj.Strains = zeros([size(obj.Points) obj.RefFrame]) ;
                    end
        end
        
        
        function obj = modify(obj,hd)
        end
        
        
        function obj = computeDisplacements(obj,hd)
            if ~obj.compDisp ; return ; end
            script = ['navDIC_',obj.displMethod] ;
            eval(['[obj,hd] = ',script,'(obj,hd) ;']) ;
        end
        
        
        function obj = computeStrains(obj,hd)
            if ~obj.compStrains ; return ; end
            script = ['navDIC_',obj.strainMethod] ;
            eval(['[obj,hd] = ',script,'(obj,hd) ;']) ;
        end
        
        
        function updateSeedPreview(obj,ax)
        end
        
    end

    
end