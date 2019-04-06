classdef navDICSlicingTool < navDICCameraPreview
    
    properties
        % CAMERA PREVIEW
            % ImLine object representing the slice
                SliceLine = [] ;
                SlicePlot = [] ;
                SliceCallBack = [] ;
        % SLICE PREVIEW
            SlicePrev = [] ;
            SliceAxes = [] ;
            SliceImg = [] ;
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
                            posLine = [nJ/2 1 ; nJ/2 nI] ;
                        else
                            posLine = [1 nI/2 ; nJ nI/2] ;
                        end
                        prev.SliceLine = imline(prev.AxesImg,posLine(:,1),posLine(:,2)) ;
                    % Constraint on the position of the line
                        fcn = makeConstrainToRectFcn('imline',[1 nJ],[1 nI]) ;
                        setPositionConstraintFcn(prev.SliceLine,fcn) ;
                    % Line Position Callback
                        prev.SliceCallBack = addNewPositionCallback(prev.SliceLine,@(pos)disp(pos)) ;
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
                    removeNewPositionCallback(prev.SliceLine,prev.SliceCallBack) ;
                    prev.SliceCallBack = addNewPositionCallback(prev.SliceLine,sliceFunction) ;
            end
            
        % SLICING FUNCTION
            function slice(prev,hd)
                % Get the Images
                    IMG = hd.Images{prev.CameraID} ;
                    [nI,nJ,nColors,nFrames] = size(IMG) ;
                    nPix = nI*nJ ;
                    if nColors>1
                        IMG = sum(IMG/nColors,3,'native') ;
                    end
                % Get the line position
                    pos = getPosition(prev.SliceLine) ;
                    L = sqrt(sum(diff(pos,1,1).^2)) ;
                % Discretize the slicing line
                    nPts = ceil(1*L) ;
                    jj = linspace(pos(1,1),pos(2,1),nPts) ;
                    ii = linspace(pos(1,2),pos(2,2),nPts) ;
                % Linear Interpolation between the four neightbors
                    % Indices corresponding to the neightboring points
                        indN1 = floor(ii) + (floor(jj)-1)*nI ;
                        indN2 = min(indN1+nI,nPix) ;
                        indN3 = min(indN2+1,nPix) ;
                        indN4 = min(indN1+1,nPix) ;
                    % Corresponding "images"
                        matFr = ones(nPts,1)*((0:nFrames-1)*nPix) ;
                        img1 = double(IMG(indN1(:)*ones(1,nFrames)+matFr)) ;
                        img2 = double(IMG(indN2(:)*ones(1,nFrames)+matFr)) ;
                        img3 = double(IMG(indN3(:)*ones(1,nFrames)+matFr)) ;
                        img4 = double(IMG(indN4(:)*ones(1,nFrames)+matFr)) ;
                    % Interpolation
                        di = ii-floor(ii) ; dj = jj-floor(jj) ;
                        DATA = bsxfun(@times,img1,(1-dj(:)).*(1-di(:))) ...
                               + bsxfun(@times,img2,(dj(:)).*(1-di(:))) ...
                               + bsxfun(@times,img3,(dj(:)).*(di(:))) ...
                               + bsxfun(@times,img4,(1-dj(:)).*(di(:))) ...
                               ;
                       %cla(prev.SliceAxes) ; plot(prev.SliceAxes,DATA(:,hd.CurrentFrame)) ;
               % Display Slice
                    xdata = 1:nFrames ;
                    ydata = 1:nPts ;
                    cdata = DATA ;
                    if isempty(prev.SliceImg)
                    else
                        prev.SliceImg.XData = xdata ;
                        prev.SliceImg.YData = ydata ;
                        prev.SliceImg.CData = cdata ;
                    end
               % Plot on the image
                   if 1
                        tangent = diff(pos,1,1)/L ; tangent = tangent/norm(tangent) ;
                        normal = [-tangent(2) tangent(1)] ;
                        data = DATA(:,hd.CurrentFrame) ;
                        data = (data - mean(data))/range(data)*min(nI,nJ)/2 ;
                        prev.SlicePlot.XData = jj(:)+normal(1)*data(:) ;
                        prev.SlicePlot.YData = ii(:)+normal(2)*data(:) ;
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