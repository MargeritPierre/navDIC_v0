classdef navDICCameraPreview < navDICPreview
    
    properties
        % Camera Infos
            CameraName = [] ;
            CameraID = [] ;
        % Axes where the camera image is displayed
            AxesImg = [] ;
        % Handles to the Image Display
            Img = [] ;
            imgMenu = gobjects(0) ;
        % Options
            ShowImage = true ;
            BlackAndWhiteImage = true ;
            Normalize = false ;
            Invert = false ;
            ShowTimeStamp = 'none' ;
            TimeStampText = gobjects(0) ;
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
                    prev.AxesImg = axes('nextplot','add','outerposition',[0 0 1 1]) ;
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
                        sz = size(hd.Images{ID}{hd.CurrentFrame}) ; 
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
                % Initialize the timestamp text
                    xPos = prev.AxesImg.XLim(1)+0.01*range(prev.AxesImg.XLim) ;
                    yPos = prev.AxesImg.YLim(1)+0.01*range(prev.AxesImg.YLim) ;
                    prev.TimeStampText = text(prev.AxesImg,xPos,yPos ...
                                                ,'' ...
                                                ,'horizontalalignment','left' ...
                                                ,'verticalalignment','top' ...
                                                ,'color','w' ...
                                                ,'FontSize',20 ...
                                                ,'Interpreter','tex' ...
                                                ,'tag','timestamp' ...
                                             ) ;
                % Image processing menu
                    submenus = gobjects(0) ;
                    prev.imgMenu = uimenu(prev.fig,'Label','Image') ;
                    % Show image
                        submenus(end+1) = uimenu(prev.imgMenu,'Label','Show','checked',prev.ShowImage) ;
                    % Coloring
                        mColor = uimenu(prev.imgMenu,'Label','Color') ;
                            submenus(end+1) = uimenu(mColor,'Label','Black&White','checked',prev.BlackAndWhiteImage) ;
                            submenus(end+1) = uimenu(mColor,'Label','RGB or Colormap','checked',~prev.BlackAndWhiteImage) ;
                    % Process
                        submenus(end+1) = uimenu(prev.imgMenu,'Label','Normalize') ;
                        submenus(end+1) = uimenu(prev.imgMenu,'Label','Invert') ;
                    % TimeStamp
                        mStamp = uimenu(prev.imgMenu,'Label','Time Stamp') ;
                            submenus(end+1) = uimenu(mStamp,'Label','none','checked','on') ;
                            submenus(end+1) = uimenu(mStamp,'Label','abs') ;
                            submenus(end+1) = uimenu(mStamp,'Label','rel') ;
                        set(submenus,'callback',@(src,evt)prev.updateMenus(src)) ;
                % Constrain the aspect ratio
                    %prev.fig.SizeChangedFcn = @(fig,evt)navDICCameraPreview.fixAspectRatio(fig,ratios) ;
                % Update the preview
                    %prev = updatePreview(prev,hd) ;
            end
            
        % UPDATE
            function prev = updateMenus(prev,submenu)
                % toggle checking
                    if submenu.Parent==prev.imgMenu % one option
                        submenu.Checked = ~submenu.Checked ;
                    else % list of options
                        % Uncheck all subMenuItems
                            set(submenu.Parent.Children,'checked','off')
                        % Check the selected item
                            submenu.Checked = 'on' ;
                    end
                % Set preview properties
                    prev.ShowImage = get(findobj(prev.imgMenu,'Label','Show'),'checked') ;
                    prev.BlackAndWhiteImage = get(findobj(prev.imgMenu,'Label','Black&White'),'checked') ;
                    prev.ShowTimeStamp = get(findobj(get(findobj(prev.imgMenu,'Label','Time Stamp'),'children'),'checked','on'),'Label') ;
                    prev.Normalize = get(findobj(prev.imgMenu,'Label','Normalize'),'checked') ;
                    prev.Invert = get(findobj(prev.imgMenu,'Label','Invert'),'checked') ;
            end
            
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
                    prev.Img.Visible = prev.ShowImage ; 
                    if prev.BlackAndWhiteImage && size(img,3) == 1 && isinteger(img) 
                        prev.Img.CData = repmat(img,[1 1 3]) ;
                    else
                        prev.Img.CData = img ;
                    end
                    if prev.Normalize
                        ranClass = getrangefromclass(prev.Img.CData) ;
                        ranImg = [min(prev.Img.CData(:)) max(prev.Img.CData(:))] ;
                        prev.Img.CData = (prev.Img.CData-ranImg(1))*(ranClass(2)/diff(ranImg))+ranClass(1) ;
                    end
                    if prev.Invert
                        ranClass = getrangefromclass(prev.Img.CData) ;
                        prev.Img.CData = ranClass(2)-prev.Img.CData ;
                    end
                % Actualize timestamp if needed
                    switch prev.ShowTimeStamp
                        case 'abs'
                            prev.TimeStampText.String = string(datetime(hd.TimeLine(hd.CurrentFrame,:))) ;
                        case 'rel'
                            prev.TimeStampText.String = string(datetime(hd.TimeLine(hd.CurrentFrame,:))-datetime(hd.TimeLine(1,:))) ;
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
        
        methods
            function set.ShowTimeStamp(prev,val)
                val = lower(val) ;
                switch val
                    case 'none'
                        prev.TimeStampText.String = '' ;
                    case 'abs'
                    case 'rel'
                    otherwise
                        error('unknown value: must be ''none'', ''abs'' or ''rel''.') ;
                end
                prev.ShowTimeStamp = val ;
            end
        end
    
end