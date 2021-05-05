classdef navDICSlicingTool < navDICCameraPreview
    
    properties
        % CAMERA PREVIEW
            % ImLine object representing the slice
                SliceLines = imline.empty ;
                SliceContour = [] ;
                SlicePlot = [] ;
                SliceCallBacks = [] ;
        % SLICE PREVIEW
            SlicePrev = [] ;
            SliceAxes = [] ;
            SliceImg = [] ;
            ModeMenus = [] ;
    end
    
    methods
        
        % CONSTRUCTOR        
            function prev = navDICSlicingTool(hd,varargin)
                % CAMERA PREVIEW
                    % Superclass constructor call (Camera Preview)
                        prev = prev@navDICCameraPreview(hd,varargin{:}) ;
                        prev.fig.Name = 'navDIC Slicing Preview' ;
                    % Lock the axes limits
                        [nI,nJ] = size(prev.Img.CData(:,:,1)) ;
                        prev.AxesImg.XLimMode = 'manual' ;
                        prev.AxesImg.XLim = [1 nJ] ;
                        prev.AxesImg.YLimMode = 'manual' ;
                        prev.AxesImg.YLim = [1 nI] ;
                    % Slice Plot
                        prev.SlicePlot = plot(prev.AxesImg,NaN,NaN,'r','linewidth',1) ;
                    % Add the slice line
                        if nI>nJ
                            posLine = [nJ/2 2 ; nJ/2 nI-1] ;
                        else
                            posLine = [2 nI/2 ; nJ-1 nI/2] ;
                        end
                        prev.SliceContour = plot(prev.AxesImg,posLine([1 2 2 1 1],1),posLine([1 1 2 2 1],2),'b') ;
                        prev.SliceLines(1) = imline(prev.AxesImg,posLine(:,1),posLine(:,2)) ;
                        prev.SliceLines(2) = imline(prev.AxesImg,posLine(:,1),posLine(:,2)) ;
                    % Constraint on the position of the line and dummy
                    % callback
                        fcn = makeConstrainToRectFcn('imline',[2 nJ-1],[2 nI-1]) ;
                        for s = 1:numel(prev.SliceLines)
                            setPositionConstraintFcn(prev.SliceLines(s),fcn) ;
                            if s==1
                                prev.SliceCallBacks = addNewPositionCallback(prev.SliceLines(s),@(pos)disp(pos)) ;
                            else
                                prev.SliceCallBacks(s) = addNewPositionCallback(prev.SliceLines(s),@(pos)disp(pos)) ;
                            end
                        end
                    % Add menus
                        % Mode menu
                            MainModeMenu = uimenu(prev.fig,'Label','Mode') ;
                                prev.ModeMenus = uimenu(MainModeMenu,'Label','Single Line','Checked','on') ;
                                prev.ModeMenus(end+1) = uimenu(MainModeMenu,'Label','Ruled Surface') ;
                                set(prev.ModeMenus,'Callback',@(src,evt)prev.MenusCallbackFcn(src,evt)) ;
                % SLICE PREVIEW
                    prev.SlicePrev = navDICPreview(hd,varargin{:}) ;
                        prev.SlicePrev.fig.Name = 'navDIC Slice Image' ;
                    prev.SliceAxes = axes() ;
                        prev.SliceAxes.YDir  = 'reverse' ;
                        axis tight
                        colormap(prev.SliceAxes,'gray') ;
                    prev.SliceImg = imagesc(rand(10)*NaN,'tag','sliceImg') ;
                % Bring the first figure to front
                    figure(prev.fig) ;
                % Set the closing function
                    prev.fig.CloseRequestFcn = @(src,evt)closePreview(prev) ;
                    prev.SlicePrev.fig.CloseRequestFcn = @(src,evt)closePreview(prev) ;
                % Update the preview
                    prev = updatePreview(prev,hd) ;
            end
            
        % UPDATE
            function prev = updatePreview(prev,hd)
                % Superclass updating function
                    prev = updatePreview@navDICCameraPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
                % Update the slice function
                    sliceFunction = @(pos)slice(prev,hd) ;
                % Execute it
                    sliceFunction() ;
                % Set it as the slice line position callback
                    for s = 1:numel(prev.SliceLines)
                        removeNewPositionCallback(prev.SliceLines(s),prev.SliceCallBacks(s)) ;
                        prev.SliceCallBacks(s) = addNewPositionCallback(prev.SliceLines(s),sliceFunction) ;
                    end
            end
            
        % SLICING FUNCTION
            function slice(prev,hd)
                disp(['slice! ',num2str(rand(1))])
                %return ;
                % Get the lines position
                    X = [] ; Y = [] ; Ls = [] ;
                    for s = 1:numel(prev.SliceLines)
                        pos = getPosition(prev.SliceLines(s)) ;
                        X(:,s) = pos(:,1) ; Y(:,s) = pos(:,2) ;
                    end
                    meanPos = [mean(X,2) mean(Y,2)] ;
                    % Slice Length
                        Ls = ceil(max(sqrt(max(diff(X,1,1).^2 + diff(Y,1,1).^2)),1)) ;
                    % Rule length
                        Lr = ceil(max(sqrt(max(diff(X,1,2).^2 + diff(Y,1,2).^2)),1)) ;
                % Depending on the mode
                    mode = get(findobj(prev.ModeMenus,'Checked','on'),'Label') ;
                    switch mode
                        case 'Single Line'
                            % Avoid looping
                            if isequal(getPosition(prev.SliceLines(1)),getPosition(prev.SliceLines(2))) ; return ; end
                            % Superimpose the two lines
                            setPosition(prev.SliceLines(1),getPosition(prev.SliceLines(2))) ;
                        case 'Ruled Surface'
                    end
                % Slice Contour
                    prev.SliceContour.XData = [X(1,1) X(1,2) X(2,2) X(2,1) X(1,1)] ;
                    prev.SliceContour.YData = [Y(1,1) Y(1,2) Y(2,2) Y(2,1) Y(1,1)] ;
                % Discretize the slicing region (Q4 FE)
                    [ee,nn] = meshgrid(linspace(0,1,Ls),linspace(0,1,Lr)) ;
                    jj = X(1,1).*(1-ee).*(1-nn) + X(2,1).*(ee).*(1-nn) + X(1,2).*(1-ee).*(nn) + X(2,2).*(ee).*(nn) ;
                    ii = Y(1,1).*(1-ee).*(1-nn) + Y(2,1).*(ee).*(1-nn) + Y(1,2).*(1-ee).*(nn) + Y(2,2).*(ee).*(nn) ;
                    %delete(findobj(prev.AxesImg,'tag','markers')) ; %pl = plot(prev.AxesImg,jj(:),ii(:),'.r','tag','markers','hittest','off') ;
                % Get the Image Infos
                    [nI,nJ,nColors] = size(hd.Images{prev.CameraID}{end}) ;
                    nFrames = numel(hd.Images{prev.CameraID}) ;
                % Real-valued pixel coord: Linear Interpolation between the four neightbors
                    % Indices corresponding to the neightboring pixels
                    % (frame space)
                        indN1 = floor(ii) + (floor(jj)-1)*nI ;
                        indN2 = indN1+nI ;
                        indN3 = indN2+1 ;
                        indN4 = indN1+1 ;
                    % Local coordinates
                        di = ii-floor(ii) ; dj = jj-floor(jj) ;
                % GET THE SLICE
                    if 0 % Small data version
                    % Transform to space-time
                        IMG = cat(4,hd.Images{prev.CameraID}{:}) ;
                        if nColors>1 ; IMG = sum(IMG/nColors,3,'native') ; end
                        IMG = reshape(IMG,[],nFrames) ; % critical step (memory needed)
                        % Interpolation
                            %SLICE = T*double(IMG) ;
                            fr = reshape((0:nFrames-1)*nI*nJ,[1 1 nFrames]) ;
                            im1 = IMG(indN1+fr) ;
                            im2 = IMG(indN2+fr) ;
                            im3 = IMG(indN3+fr) ;
                            im4 = IMG(indN4+fr) ;
                            SLICE = double(im1).*(1-di).*(1-dj) ...
                                    + double(im2).*(dj).*(1-di) ...
                                    + double(im3).*(dj).*(di) ...
                                    + double(im4).*(1-dj).*(di) ...
                                    ;
                            SLICE = squeeze(mean(SLICE,1)) ;
                    elseif 1 % Large data version
                    % Build the transfer matrix
                        iit = repmat(1:numel(indN1),[1 4])' ;
                        jjt = [indN1(:);indN2(:);indN3(:);indN4(:)] ;
                        di = di(:) ; dj = dj(:) ;
                        vvt = [(1-di).*(1-dj);(dj).*(1-di);(dj).*(di);(1-dj).*(di)] ;
                        T = sparse(iit(:),jjt(:),vvt(:),Lr*Ls,nI*nJ) ;
                    % Extract the slice
                        SLICE = cellfun(@(img)T*reshape(double(img),[nI*nJ 1 1 nColors]),hd.Images{prev.CameraID},'UniformOutput',false) ; 
                        SLICE = reshape(cat(3,SLICE{:}),[Lr Ls nFrames nColors]) ;
                        SLICE = reshape(mean(SLICE,1),[Ls nFrames nColors]) ;
                    end
               % Display Slice
                    xdata = 1:nFrames ;
                    ydata = 1:Ls ;
                    cdata = SLICE ;
                    if isempty(prev.SliceImg)
                    else
                        prev.SliceImg.XData = xdata ;
                        prev.SliceImg.YData = ydata ;
                        prev.SliceImg.CData = cdata ;
                    end
               % Plot on the image
                   if 1
                        tangent = diff(meanPos,1,1)/Ls ; tangent = tangent/norm(tangent) ;
                        normal = [-tangent(2) tangent(1)] ;
                        data = SLICE(:,hd.CurrentFrame) ;
                        data = (data - mean(data))/range(data)*min(nI,nJ)/2 ;
                        xx = linspace(meanPos(1,1),meanPos(2,1),Ls) ;
                        yy = linspace(meanPos(1,2),meanPos(2,2),Ls) ;
                        prev.SlicePlot.XData = xx(:)+normal(1)*data(:) ;
                        prev.SlicePlot.YData = yy(:)+normal(2)*data(:) ;
                   end
            end
            
        % MENU UPDATING
            function MenusCallbackFcn(prev,src,evt)
                switch src.Parent.Label
                    case 'Mode'
                        set(src.Parent.Children,'Checked','off') ;
                        src.Checked = 'on' ;
                        setPosition(prev.SliceLines(1),getPosition(prev.SliceLines(2))) ;
                end
            end
        
        % DESTRUCTOR
            function closePreview(obj)
                figsToClose = [obj.fig obj.SlicePrev.fig] ;
                set(figsToClose,'CloseRequestFcn','closereq') ;
                close(figsToClose) ;
            end
    end
    
            
    % OTHER FUNCTIONS
    
        methods(Static)
                
        end
    
end