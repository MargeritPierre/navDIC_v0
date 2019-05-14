classdef navDICCameraPreview < navDICPreview
    
    properties
        % Camera Infos
            CameraName = [] ;
            CameraID = [] ;
        % Axes where the image is displayed
            AxesImg = [] ;
        % Handles to the Image Display
            Img = [] ;
    end
    
    methods
        
        % CONSTRUCTOR        
            function prev = navDICCameraPreview(hd,varargin)
                % Superclass constructor call
                    prev = prev@navDICPreview(hd,varargin{:}) ;
                    prev.fig.Name = 'navDIC Camera Preview' ;
                % Choose a camera to preview
                    if ~isempty(varargin)
                        ID = varargin{1} ;
                        valid = ID~=-1 ;
                    else
                        [ID,valid] = selectCameras(hd,'single') ;
                    end
                    if ~valid 
                        close(prev.fig)
                        prev.isValid = false ;
                        return ;
                    end
                % Retrieve Camera Infos
                    prev.CameraName = hd.Cameras(ID).Name ;
                    prev.CameraID = ID ;
                    prev.fig.Name = [prev.fig.Name,' : ',prev.CameraName] ;
                % Set the axes
                    prev.AxesImg = axes('position',[0 0 1 1]) ;
                        axis tight
                        axis off
                        axis equal
                        prev.AxesImg.YDir = 'reverse' ;
                % Set the figure at the right size
                    prev.fig.Units = 'pixels' ;
                    posFig = prev.fig.Position(3:4) ;
                    roiCam = hd.Cameras(ID).VidObj.ROIPosition ;
                    resCam = roiCam(3:4) ;
                    ratios = resCam./posFig ;
                    ratios = ratios./max(ratios) ;
                    newSizeFig = posFig.*ratios ;
                    prev.fig.Position = [prev.fig.Position(1:2)+(posFig-newSizeFig)/2 newSizeFig] ;
                % Create an empty image
                    prev.Img = image(ones([fliplr(resCam) 3])*.5) ;
                % Contrain the aspect ratio
                    %prev.fig.SizeChangedFcn = @(fig,evt)navDICCameraPreview.fixAspectRatio(fig,ratios) ;
            end
            
        % UPDATE
            function prev = updatePreview(prev,hd)
                % Superclass updating function
                    prev = updatePreview@navDICPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
                % Try to get the last acquired image on camera
                    img = [] ;
                    try
                        img = hd.Images{hd.CurrentFrame} ;
                        img = img{prev.CameraID} ;
                    end
                    if isempty(img) 
                        prev.Img.CData = 0.5 + 0*prev.Img.CData ;
                        return ;
                    end
                % Actualize preview image
                    nBands = size(img,3) ;
                    if nBands == 1
                        while iscell(img)
                            img = img{1};
                        end
                        prev.Img.CData = repmat(img,[1 1 3]) ;
                    else

                        prev.Img.CData = img ;
                    end
            end
        
        % DESTRUCTOR
            function closePreview(obj)
            end
    end
    
            
    % OTHER FUNCTIONS
    
        methods(Static)

            % FIXED ASPECT RATIO
                function fixAspectRatio(fig,fixedRatio)
                    pause(0.02) ;
                    fig.Units = 'normalized' ;
                    posFig = fig.Position ;
                    sizeFig = posFig(4).*fixedRatio ;
                    fig.Position = [posFig(1:2)+(posFig(3:4)-sizeFig)/2 sizeFig] ;
                end
                
        end
    
end