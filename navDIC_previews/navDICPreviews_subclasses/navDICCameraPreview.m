classdef navDICCameraPreview < navDICPreview
    
    properties
        % Camera Infos
            CameraName = [] ;
            CameraID = [] ;
        % Axes where the camera image is displayed
            AxesImg = [] ;
        % Handles to the Image Display
            Img = [] ;
        % Options
            BlackAndWhiteImage = true ;
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
                    prev.AxesImg = axes('outerposition',[0 0 1 1]) ;
                        prev.AxesImg.YDir = 'reverse' ;
                        prev.AxesImg.XTick = [] ;
                        prev.AxesImg.YTick = [] ;
                        prev.AxesImg.LooseInset = [1 1 1 1]*0 ;
                        prev.AxesImg.Clipping = 'off' ;
                        axis tight
                        axis off
                        axis equal
                % Set the figure at the right size
                    prev.fig.Units = 'pixels' ;
                    posFig = prev.fig.Position(3:4) ;
                    if ~isempty(hd.Images) && ~isempty(hd.Images{ID})
                        sz = size(hd.Images{ID}{end}) ; 
                        resCam = sz([2,1]) ;
                    else
                        roiCam = hd.Cameras(ID).VidObj.ROIPosition ; 
                        resCam = roiCam(3:4) ;
                    end
                    ratios = resCam./posFig ;
                    ratios = ratios./max(ratios) ;
                    newSizeFig = posFig.*ratios ;
                    prev.fig.Position = [prev.fig.Position(1:2)+(posFig-newSizeFig)/2 newSizeFig] ;
                % Create an empty image as surface (for caxis scaling)
                    img = ones([fliplr(resCam) 3])*.5 ;
                    prev.Img = imagesc(img) ;
                    %prev.Img = surf(JJ-0.5,II-0.5,JJ*0,img,'facecolor','flat','edgecolor','none') ;
                % Contrain the aspect ratio
                    %prev.fig.SizeChangedFcn = @(fig,evt)navDICCameraPreview.fixAspectRatio(fig,ratios) ;
                % Update the preview
                    %prev = updatePreview(prev,hd) ;
            end
            
        % UPDATE
            function prev = updatePreview(prev,hd)
                % Superclass updating function
                    prev = updatePreview@navDICPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
                % Try to get the last acquired image on camera
                    img = [] ;
                    try
                        img = hd.Images{prev.CameraID}{hd.CurrentFrame} ;
                    end
                    if isempty(img) 
                        %prev.Img.CData = 0.5 + 0*prev.Img.CData ;
                        return ;
                    end
                % Actualize preview image
                    % Process
                        %img = single(img) ;
                        %img = img*(max(getrangefromclass(img(:)))/range(img(:))) ;
                    if prev.BlackAndWhiteImage && size(img,3) == 1 && isinteger(img) 
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