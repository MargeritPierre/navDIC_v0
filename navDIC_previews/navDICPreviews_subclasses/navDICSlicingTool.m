classdef navDICSlicingTool < navDICCameraPreview
    
    properties
        % CAMERA PREVIEW
            % ImLine object representing the slice
                SliceLines = imline.empty ;
                SliceContour = [] ;
                SlicePlot = [] ;
                SliceCallBacks = [] ;
            % Menus
                ModeMenus = [] ;
                PresetMenus = [] ;
                ViewMenus = [] ;
        % SLICE PREVIEW
            SlicePrev = [] ;
            SliceAxes = [] ;
            SliceImg = [] ;
            FrameMarker = [] ;
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
                        % Preset menus
                            MainPresetMenu = uimenu(prev.fig,'Label','Preset') ;
                                prev.PresetMenus = uimenu(MainPresetMenu,'Label','Horizontal','Checked','on') ;
                                prev.PresetMenus(end+1) = uimenu(MainPresetMenu,'Label','Vertical') ;
                                set(prev.PresetMenus,'Callback',@(src,evt)prev.MenusCallbackFcn(src,evt)) ;
                % SLICE PREVIEW
                    prev.SlicePrev = navDICPreview(hd,varargin{:}) ;
                        prev.SlicePrev.fig.Name = 'navDIC Slice Image' ;
                    prev.SliceAxes = axes('nextplot','add') ;
                        prev.SliceAxes.YDir  = 'reverse' ;
                        axis tight
                        colormap(prev.SliceAxes,'gray') ;
                    prev.SliceImg = imagesc(rand(10)*NaN,'tag','sliceImg') ;
                    prev.FrameMarker = plot(prev.SliceAxes,NaN,NaN,'r','linewidth',1) ;
                % View menus
                    MainViewMenu = uimenu(prev.SlicePrev.fig,'Label','View') ;
                        prev.ViewMenus = uimenu(MainViewMenu,'Label','Fixed','Checked','on') ;
                        prev.ViewMenus(end+1) = uimenu(MainViewMenu,'Label','Left') ;
                        prev.ViewMenus(end+1) = uimenu(MainViewMenu,'Label','Center') ;
                        prev.ViewMenus(end+1) = uimenu(MainViewMenu,'Label','Right') ;
                        set(prev.ViewMenus,'Callback',@(src,evt)prev.MenusCallbackFcn(src,evt)) ;
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
                % Set it as the slice line position callback
                    for s = 1:numel(prev.SliceLines)
                        removeNewPositionCallback(prev.SliceLines(s),prev.SliceCallBacks(s)) ;
                        prev.SliceCallBacks(s) = addNewPositionCallback(prev.SliceLines(s),sliceFunction) ;
                    end
                % IF the number of frames has changed, change the slice image
                    if size(prev.SliceImg.CData,2)~=hd.nFrames
                        sliceFunction([]) ;
                    else % Update the slice plot only
                        prev.updatePlot(hd) ;
                    end
               
            end
            
        % SLICING FUNCTION
            function slice(prev,hd)
                disp(['slice! ',num2str(clock()*[0;0;0;3600;60;1])])
                dataSize = numel(hd.Images{prev.CameraID}{1}(:))*numel(hd.Images{prev.CameraID}) ;
                smallData = dataSize<1e10 ;
                useTransferMat = false ;
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
                            % Superimpose the two lines
                            if ~isequal(getPosition(prev.SliceLines(1)),getPosition(prev.SliceLines(2)))
                                setPosition(prev.SliceLines(1),getPosition(prev.SliceLines(2))) ; 
                                return ; % Avoid recursive looping
                            end
                        case 'Ruled Surface'
                    end
                % Slice Contour
                    prev.SliceContour.XData = [X(1,1) X(1,2) X(2,2) X(2,1) X(1,1)] ;
                    prev.SliceContour.YData = [Y(1,1) Y(1,2) Y(2,2) Y(2,1) Y(1,1)] ;
                % Discretize the slicing region (Q4 FE)
                    [ee,nn] = meshgrid(linspace(0,1,Ls),linspace(0,1,Lr)) ;
                    ee = ee(:) ; nn = nn(:) ;
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
                % PREPARE IMAGE DATA
                    if smallData % Small data version
                    % Transform to space-time array
                        IMG = {cat(4,hd.Images{prev.CameraID}{:})} ; % memory critical ! {[nI nJ nColors nFrames]}
                    else % Large data version
                    % Keep the cell array format
                        IMG = hd.Images{prev.CameraID} ;
                    end
                % Convert to double
                    for fr = 1:numel(IMG) ; IMG{fr} = reshape(double(IMG{fr}),nI*nJ,[]) ; end
                % GET THE SLICE
                    SLICE = cell(numel(IMG),1) ;
                    if useTransferMat % Use a transfer matrix
                    % Build the transfer matrix
                        iit = repmat(1:numel(indN1),[1 4])' ;
                        jjt = [indN1(:);indN2(:);indN3(:);indN4(:)] ;
                        di = di(:) ; dj = dj(:) ;
                        vvt = [(1-di).*(1-dj);(dj).*(1-di);(dj).*(di);(1-dj).*(di)] ;
                        T = sparse(iit(:),jjt(:),vvt(:),Lr*Ls,nI*nJ) ;
                    % Project
                        for fr = 1:numel(IMG) ; SLICE{fr} = reshape(T*IMG) ; end
                    else % Use indexing
                        for fr = 1:numel(IMG)
                            SLICE{fr} = IMG{fr}(indN1,:).*(1-di).*(1-dj) ...
                                    + IMG{fr}(indN2,:).*(dj).*(1-di) ...
                                    + IMG{fr}(indN3,:).*(dj).*(di) ...
                                    + IMG{fr}(indN4,:).*(1-dj).*(di) ...
                                    ;
                        end
                    end
                % Concatenate frames (if needed)
                    SLICE = cat(3,SLICE{:}) ;
                % Compute the mean along the rule length
                    SLICE = mean(reshape(SLICE,[Lr Ls nColors nFrames]),1) ; % [1 Ls nColors nFrames]
                % Permute dimensions
                    SLICE = permute(SLICE,[2 4 3 1]) ;
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
                   updatePlot(prev,hd) ;
            end
            
        % Plot the slice at the current frame
            function updatePlot(prev,hd)                
            % Get the Image Infos
                [nI,nJ,~] = size(hd.Images{prev.CameraID}{end}) ;
            % Get the lines position
                X = [] ; Y = [] ; Ls = [] ;
                for s = 1:numel(prev.SliceLines)
                    pos = getPosition(prev.SliceLines(s)) ;
                    X(:,s) = pos(:,1) ; Y(:,s) = pos(:,2) ;
                end
                meanPos = [mean(X,2) mean(Y,2)] ;
            % Slice Length
                Ls = ceil(max(sqrt(max(diff(X,1,1).^2 + diff(Y,1,1).^2)),1)) ;
            % Line normal/tangent
                tangent = diff(meanPos,1,1)/Ls ; tangent = tangent/norm(tangent) ;
                normal = [-tangent(2) tangent(1)] ;
            % Retrieve data
                data = mean(prev.SliceImg.CData(:,hd.CurrentFrame,:),3) ;
            % Scale data
                amp = min(nI,nJ)./2.5 ;
                data = data./range(data,1).*amp ;
                data = data*2-amp ;
                %data = data - mean(data,1) ;
            % Plot position
                xx = linspace(meanPos(1,1),meanPos(2,1),Ls) ;
                yy = linspace(meanPos(1,2),meanPos(2,2),Ls) ;
            % Update the slice plot
                prev.SlicePlot.XData = xx(:)+normal(1)*data(:) ;
                prev.SlicePlot.YData = yy(:)+normal(2)*data(:) ;
            % Update the frame marker
                prev.FrameMarker.XData = [1 1]*hd.CurrentFrame ;
                prev.FrameMarker.YData = [min(prev.SliceImg.YData(:)) max(prev.SliceImg.YData(:))] ;
            % Update the view limits
                Lf = range(prev.SliceAxes.XLim)-1 ;
                xx = prev.SliceImg.XData(:) ; 
                switch prev.ViewMenus([prev.ViewMenus.Checked]).Label
                    case 'Fixed'
                        xlim = [] ;
                    case 'Left'
                        xlim = min(hd.CurrentFrame,max(xx)-Lf) + [0 Lf] ;
                    case 'Center'
                        xlim = max(min(hd.CurrentFrame,max(xx)-Lf/2),min(xx)+Lf/2) + [-Lf Lf]/2 ;
                    case 'Right'
                        xlim = max(hd.CurrentFrame,min(xx)+Lf) + [-Lf 0] ;
                end
                if ~isempty(xlim) ; prev.SliceAxes.XLim = xlim + [-1 1]*0.5 ; end
            end
            
        % MENU UPDATING
            function MenusCallbackFcn(prev,src,evt)
            % React to a change in the menus
            % Change the checked menu
                set(src.Parent.Children,'Checked','off') ;
                src.Checked = 'on' ;
            % Apply changes
                switch src.Parent.Label
                    case 'Mode'
                        setPosition(prev.SliceLines(1),getPosition(prev.SliceLines(2))) ;
                    case 'Preset'
                        [nI,nJ] = size(prev.Img.CData(:,:,1)) ;
                        switch src.Label
                            case 'Vertical'
                                setPosition(prev.SliceLines(2),[nJ/2 2 ; nJ/2 nI-1]) ;
                            case 'Horizontal'
                                setPosition(prev.SliceLines(2),[2 nI/2 ; nJ-1 nI/2]) ;
                        end
                    case 'View'
                        setPosition(prev.SliceLines(2),getPosition(prev.SliceLines(2))) ;
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