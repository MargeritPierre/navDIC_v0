function [setup,hd] = loadVideoFrames(hd,filename)
            
    % Initialize the setup
        [path,file,ext] = fileparts(filename) ;
        setup.Valid = false ;
        setup.Path = path ;
        setup.CommonName = file ;
        setup.ImagesExtension = ext ;
        
    % Initialize the nested variables
        video = [] ;
        H = [] ;
        
    % Try to load the file
        disp('------ VIDEO LOADING ------')
        disp(['   File: ',filename])
        disp(['   Loading...']) ;
        video = VideoReader(filename) ;
        
    % Display Video Informations
        get(video)
        
    % Prompt the figure to choose options
        initOptionsFigure ;
        
    % Wait for the user to close the figure
        while isvalid(H.fig) && ~H.btnOK.Value
            drawnow ;
        end
        if ~isvalid(H.fig) ; return ; end % The figure has been closed, cancel the operation
    
    % Ask the User one last time
        if H.estimSizeOnMemory>1e8 % More than 100MB of Data wil be loaded
            answer = questdlg({['The estimated data size to be loaded is ',num2str(H.estimSizeOnMemory*1e-6,'%.0f'),'MB'],...
                                'Are you sure you want to proceed ?'}...
                                ,'WARNING','OK','Cancel','Cancel'...
                                ) ;
            if strcmp(answer,'Cancel') ; return ; end
            drawnow ;
        end
        
    % Load the images
        IMG = zeros([size(H.processImg(:,:,1)) size(H.processImg,3) length(H.loadedFrames)],class(H.processImg)) ;
        wtbr = waitbar(0,'Loading Frames...') ;
        for fr = 1:length(H.loadedFrames)
            % Load the frame
                currentImg = read(video,H.loadedFrames(fr)) ;
            % Process the frame
                for p = 1:length(H.imgProcesses)
                    currentImg = H.imgProcesses{p}(currentImg) ;
                end
            % Push it on the images
                IMG(:,:,:,fr) = currentImg ;
            % Waitbar
                wtbr = waitbar(fr/length(H.loadedFrames),wtbr,['Loading Frames... (',num2str(fr),'/',num2str(length(H.loadedFrames)),')']) ;
        end
        delete(wtbr)
        
        
    % SET THE CAMERA
        Camera.Name = file ;
        Camera.CurrentState = 'ghost' ;
        Camera.Adaptator = 'folder' ;
        Camera.VidObj.ROIPosition = [0 0 flip(size(IMG(:,:,1)))] ;

    % Change the handles
        hd.Images = {} ;
        hd.Images{1} = IMG ;
        hd.InputData = [] ;
        hd.nFrames = max(size(IMG,4)) ;
        hd.CurrentFrame = 1 ;
        hd.Cameras = Camera ;

    % DEFAULT TIMELINE
        hd.TimeLine = 1/video.FrameRate*(0:hd.nFrames-1)'*[0 0 0 0 0 1] ;
            
    % Validate the setup
        setup.Valid = true ;
        close(H.fig) ;

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% UTIL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% MAIN FIGURE FUNCTION
    function initOptionsFigure
        % Handle Structure
            H = [] ;
            H.Valid = false ;
            H.nFrames = video.NumberOfFrames ;
            H.vidRes = [video.Width video.Height] ;
        % Loaded frames
            H.startFrame = 1 ;
            H.endFrame = H.nFrames ;
            H.decimFrames = 1 ;
            H.loadedFrames = H.startFrame:H.decimFrames:H.endFrame ;
        % Take a reference image
            H.currentFrame = 1 ;
            H.currentImg = read(video,H.currentFrame) ;
            H.imgProcesses = {} ;
        % Init Figure
            imgRelSize = 0.6 ; % Relative size of the image to the screen
            uiHeight = 0.07 ; % Relative vertical size of the top slider to the figure
            menuSizeLeft = 0.17 ; % Relative horizontal size of the left menu to the figure
            margin = 0.005 ;
            H.fig = figure ;
                H.fig.ToolBar = 'none' ;
                H.fig.MenuBar = 'none' ;
                H.fig.NumberTitle = 'off' ;
                H.fig.Name = 'LOAD A VIDEO: CHOOSE THE PROCESSING PARAMETERS' ;
                % Positionning
                    ratios = H.vidRes./H.fig.Position(3:4) ;
                    dimFig = imgRelSize/max(ratios)*H.vidRes ;
                H.fig.Position = [H.fig.Position(1:2)+H.fig.Position(3:4)/2 0 0] ...
                                    + [-1/2 0 1 0]*dimFig(1)*1/(1-menuSizeLeft) ...
                                    + [0 -1/2 0 1]*dimFig(2)*1/(1-uiHeight) ;
        % Axes containing the image preview
            H.axImg = axes('outerposition',[menuSizeLeft 0 1-menuSizeLeft 1-uiHeight]) ;
                H.axImg.LooseInset = [1 1 1 1]*0.005 ;
                H.axImg.YDir = 'reverse' ;
                H.axImg.XTick = [] ; H.axImg.YTick = [] ;
                axis tight
                axis equal
            H.ImgH = image(H.currentImg) ;
        % ImRect for the ROI
            H.rectROI = imrect(H.axImg,[1 1 H.vidRes]) ;
            fcn = makeConstrainToRectFcn('imrect',[1 H.vidRes(1)],[1 H.vidRes(2)]) ;
            setPositionConstraintFcn(H.rectROI,fcn) ;
            addNewPositionCallback(H.rectROI,@(pos)Global_Callback('normal')) ;
        % Slider for the frame choice
            H.sliderFrame = uicontrol(H.fig,'style','slider','units','normalized'...
                                            ,'position',[menuSizeLeft+margin 1-uiHeight+margin 1-margin*2-menuSizeLeft uiHeight-margin*2] ...
                                            ,'min',1 ...
                                            ,'max',H.nFrames ...
                                            ,'value',1 ...
                                            ,'sliderstep',[1/H.nFrames 1/10] ...
                                        ) ;
            addlistener(H.sliderFrame,'Value','PostSet',@(src,evt)Global_Callback('normal'));
        % Frame Infos
            H.frameText = uicontrol(H.fig,'style','text','units','normalized'...
                                            ,'position',[margin 1-uiHeight+margin menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string',['Current Frame: ',num2str(H.currentFrame),'/',num2str(H.nFrames)] ...
                                            ... ,'HorizontalAlignment','right' ...
                                        ) ; 
        % Start|End Frames
            H.startBtn = uicontrol(H.fig,'style','pushbutton','units','normalized'...
                                            ,'position',[margin 1-2*(uiHeight+margin) menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Set Start Frame' ...
                                            ,'callback',@(src,evt)Global_Callback('SetStart') ...
                                        ) ; 
            H.endBtn = uicontrol(H.fig,'style','pushbutton','units','normalized'...
                                            ,'position',[margin 1-3*(uiHeight+margin) menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Set End Frame' ...
                                            ,'callback',@(src,evt)Global_Callback('SetEnd') ...
                                        ) ; 
        % Frame decimation
            H.decimSlider = uicontrol(H.fig,'style','slider','units','normalized'...
                                            ,'position',[margin 1-4*(uiHeight+margin) menuSizeLeft-02*margin uiHeight-2*margin] ...
                                            ,'callback',@(src,evt)Global_Callback('SetEnd') ...
                                            ,'min',0 ...
                                            ,'max',log10(H.nFrames) ...
                                            ,'value',0 ...
                                        ) ; 
            addlistener(H.decimSlider,'Value','PostSet',@(src,evt)Global_Callback('normal'));
        % Frame selection infos
            H.frameSelectText = uicontrol(H.fig,'style','text','units','normalized'...
                                            ,'position',[margin 1-5*(uiHeight+margin) menuSizeLeft-02*margin uiHeight-2*margin] ...
                                            ,'string',['Frames Selected: ' num2str(H.startFrame) ':' num2str(H.decimFrames) ':' num2str(H.endFrame)] ...
                                        ) ; 
        % White & Black
            H.convertWB = uicontrol(H.fig,'style','checkbox','units','normalized'...
                                            ,'position',[margin 1-6*(uiHeight+margin) menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Convert to White & Black' ...
                                            ,'Value',true ...
                                            ,'callback',@(src,evt)Global_Callback('normal') ...
                                        ) ; 
        % Validation Button
            H.btnOK = uicontrol(H.fig,'style','togglebutton','units','normalized'...
                                            ,'position',[margin margin menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Load Frames' ...
                                            ,'value',false ...
                                        ) ; 
        % Execute the callback once
            Global_Callback('normal') ;
    end
    
    
% GLOBAL CALLBACK 
    function Global_Callback(order)
        % PREREQUISITES
            % Current frame
                fr = round(H.sliderFrame.Value) ;
        % IMAGE PROCESSES
            % Update Image
                if H.currentFrame~=fr
                    H.currentFrame = fr ;
                    H.currentImg = read(video,H.currentFrame) ;
                    H.frameText.String = ['Current Frame: ',num2str(H.currentFrame),'/',num2str(H.nFrames)] ;
                end
            % Display Image
                H.ImgH.CData = H.currentImg ;
            % Declare image Processes
                nColors = size(H.currentImg,3) ;
                H.imgProcesses = {} ;
                H.imgProcesses{end+1} = @(img)imcrop(img,getPosition(H.rectROI)) ; % Cropping
                if H.convertWB.Value ; H.imgProcesses{end+1} = @(img)sum(img/nColors,3,'native') ; end % To White and Black
            % Apply to the current image
                H.processImg = H.currentImg ;
                for p = 1:length(H.imgProcesses)
                    H.processImg = H.imgProcesses{p}(H.processImg) ;
                end
        % VIDEO FRAMES
            % Start|step|end
                if strcmp(order,'SetStart') && fr<=H.endFrame ; H.startFrame = fr ; end
                if strcmp(order,'SetEnd') && fr>=H.startFrame ; H.endFrame = fr ; end
                H.decimFrames = round(10^H.decimSlider.Value) ;
                loadFrames = [num2str(H.startFrame) ':' num2str(H.decimFrames) ':' num2str(H.endFrame)] ;
                H.frameSelectText.String = ['Frames Selected: ' loadFrames] ;
                H.loadedFrames = eval(loadFrames) ; %H.startFrame:H.decimFrames:H.endFrame ;
            % Evaluate the total video size on memory
                img = H.processImg ;
                infos = whos('img','variables') ;
                H.estimSizeOnMemory = infos.bytes*length(H.loadedFrames) ;
                disp(['Estimated Size On Memory: ',num2str(H.estimSizeOnMemory*1e-9,'%.3f'),'GB']) ;
    end
        
      



end