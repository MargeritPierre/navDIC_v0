function [valid,hd] = loadFrames(hd,dataType,camID)

    % Initialize the output
        valid = false ;

    % Initialization
        H = [] ; % Handle structure
        H.Valid = false ;
        switch dataType
            case 'ImageFolder'
                initImageFolder ;
            case 'Video'
                initVideo ;
        end
        if ~H.Valid ; return ; end
        %disp(H) ;
        
    % Loaded frames
        H.startFrame = 1 ;
        H.endFrame = H.nFrames ;
        H.decimFrames = 1 ;
        H.loadedFrames = H.startFrame:H.decimFrames:H.endFrame ;
        
    % Take a reference image
        H.currentFrame = 1 ;
        H.currentImg = H.loadFrame(H.currentFrame) ;
        H.imgProcesses = {} ;
        
    % Prompt the figure for pre-processing
        initOptionsFigure ;
        
    % Wait for the user to close the figure
        while isvalid(H.fig) && ~H.btnOK.Value
            drawnow ;
        end
        if ~isvalid(H.fig) ; return ; end % The figure has been closed, cancel the operation
    
    % Ask the User if the data to load is large
        if H.estimSizeOnMemory>1e8 % More than 100MB of Data may be loaded
            answer = questdlg({['The estimated data size to be loaded is ',num2str(H.estimSizeOnMemory*1e-6,'%.0f'),'MB'],...
                                'Are you sure you want to proceed ?'}...
                                ,'WARNING','OK','Cancel','Cancel'...
                                ) ;
            if strcmp(answer,'Cancel') ; return ; end
            drawnow ;
        end
        
    % Load the images
        % Each frames
            IMG = cell(length(H.loadedFrames),1) ;
            wtbr = waitbar(0,'Loading Frames...') ;
            for fr = 1:length(H.loadedFrames)
                % Load the frame
                    currentImg = H.loadFrame(H.loadedFrames(fr)) ;
                % Process the frame
                    for p = 1:length(H.imgProcesses)
                        currentImg = H.imgProcesses{p}(currentImg) ;
                    end
                % Push it on the images
                    IMG{fr} = currentImg ;
                % Waitbar
                    wtbr = waitbar(fr/length(H.loadedFrames),wtbr,['Loading Frames... (',num2str(fr),'/',num2str(length(H.loadedFrames)),')']) ;
            end
        % Global Normalization
            if H.normalizeGlobal.Value
                wtbr = waitbar(1,wtbr,'Normalization...') ; drawnow ;
                Imin = min(cellfun(@(ii)min(ii(:)),IMG)) ;
                Ifactor = double(max(getrangefromclass(IMG{end})))/double(max(cellfun(@(ii)max(ii(:)),IMG))) ;
                for fr = 1:numel(IMG)
                    IMG{fr} = (IMG{fr}-Imin)*Ifactor ;
                end
            end
        delete(wtbr)
        
        
    % CHANGE THE HANDLES IF A NEW CAMERA HAS BEEN ASKED FOR
        if camID == length(hd.Cameras)+1
            % New Camera
                Camera = [] ;
                Camera.Name = H.CamName ;
                Camera.CurrentState = 'ghost' ;
                Camera.Adaptator = 'folder' ;
                Camera.VidObj.ROIPosition = [0 0 flip(size(IMG{end},[1 2]))] ;
                if camID==1
                    hd.Cameras = Camera ; % Initialize the camera list
                else
                    hd.Cameras(camID) = Camera ;
                end
            % New Default Timeline
                hd.TimeLine = H.FrameRate*(0:hd.nFrames-1)'*[0 0 0 0 0 1] ;
            % New number of frames
                hd.nFrames = max(numel(IMG),hd.nFrames) ;
        end

    % Add the images
        hd.Images{camID} = IMG ;
            
    % Validate the setup
        valid = true ;
        close(H.fig) ;
        
        
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% HANDLES INITIALIZATION FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% INITIALIZATION FUNCTIONS
    function initImageFolder
        H.Valid = false ;
        % IMAGE SOURCE
            % Choose the image folder
                [path] = uigetdir(hd.WorkDir.Path,'SELECT AN IMAGE FOLDER') ;
                if path==0 ; return ; end
                disp(newline)
                disp(['LOADING FRAMES FROM ',path])
            % Build File List
                supported = imformats ;
                imgExt = [supported.ext] ; % Image files currently supported by Matlab
                files = dir('*.anImpossibleExtension') ; % Initialization
                for i = 1:length(imgExt)
                    f = dir([path,'/*.',imgExt{i}]) ;
                    files(end+(1:length(f))) = f ;
                end
            % Keep file names only
                fileNames = {files.name} ;
                if isempty(fileNames)
                    warning(['No Valid Image Files Found in',path])
                    return
                end
            % Camera Name
                camName = strsplit(path,{'/','\'}) ;
                camName = camName{end} ;
                disp(['   CameraName: ',camName])
            % Get the common name and extension
                str = strsplit(fileNames{1},'_') ;
                if length(str)==1
                    % Empty common name ?
                        if strcmp(str,fileNames{1})
                            commonName = '' ;
                        else
                            commonName = str{1} ;
                        end
                else
                    commonName = [strjoin(str(1:end-1),'_'),'_'] ;
                end
                [~,~,ext] = fileparts(str{end}) ;
                disp(['   CommonName: ',commonName])
                disp(['   Type: ',ext])
            % Get image ids
                idSTR = {} ;
                idNUM = [] ;
                for i = 1:length(fileNames)
                    idSTR{i} = fileNames{i}(length(commonName)+1:end-length(ext)) ;
                    if ~isempty(str2double(idSTR{i}))
                        idNUM(i) = str2double(idSTR{i}) ;
                    else
                        idNUM(i) = NaN ;
                    end
                end
            % Sort images by name
                [idNUM,ind] = sort(idNUM(~isnan(idNUM))) ;
                fileNames = fileNames(ind) ;
                idSTR = idSTR(ind(~isnan(idNUM))) ;
                nFrames = length(idSTR) ;
                disp(['   Frames: [',num2str(min(idNUM)),'->',num2str(max(idNUM)),'] (',num2str(nFrames),')'])
        % IMAGE LOADING AND PROCESSING
            % Load function
                nFrames = numel(fileNames) ;
                loadFrame = @(id) imread([path,filesep,fileNames{id}]) ;
                %loadFrame = @(id) imread([path,filesep,commonName,idSTR{id},ext]) ;
                imData = loadFrame(1) ;
            % Get Infos
                [nI,nJ,nColors] = size(imData) ;
                dataType = class(imData) ;
                disp(['   Class: ',dataType])
                disp(['   Resolution: ',num2str(nJ),'x',num2str(nI)])
                disp(['   Colors: ',num2str(nColors)])
        % HANDLE STRUCTURE
            H.CamName = camName ;
            H.nFrames = nFrames ;
            H.loadFrame = loadFrame ;
            H.vidRes = [nJ nI] ;
            H.nColors = nColors ;
            H.FrameRate = 1 ;
        % Validation
            H.Valid = true ;
    end

    function initVideo
        H.Valid = false ;
        % VIDEO FILE
            % Choose the file
                supported  = getFilterSpec(VideoReader.getFileFormats()) ;
                [file,path] = uigetfile(supported,'SELECT A VIDEO FILE',[hd.WorkDir.Path,'\*']) ;
                if path==0 ; return ; end
                filename = [path file] ;
            % Try to load the file
                disp(newline) ;
                disp(['LOADING FRAMES FROM ',filename])
                video = VideoReader(filename) ;
            % Display Video Informations
                get(video)
        % LOADING FRAME FUNCTION
            % Load function
                loadFrame = @(id) read(video,id) ;
                imData = loadFrame(1) ;
            % Get Infos
                [nI,nJ,nColors] = size(imData) ;
                dataType = class(imData) ;
                disp(['   Class: ',dataType])
                disp(['   Resolution: ',num2str(nJ),'x',num2str(nI)])
                disp(['   Colors: ',num2str(nColors)])
        % HANDLE STRUCTURE
            [~,H.CamName,~] = fileparts(filename) ;
            H.nFrames = video.NumberOfFrames ;
            H.loadFrame = loadFrame ;
            H.vidRes = [nJ nI] ;
            H.nColors = nColors ;
            H.FrameRate = video.FrameRate ;
        % Validation
            H.Valid = true ;
    end
    
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% FIGURE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
% GLOBAL CALLBACK 
    function Global_Callback(order)
        % PREREQUISITES
            % Current frame
                fr = round(H.sliderFrame.Value) ;
        % IMAGE PROCESSES
            H.imgProcesses = {} ;
            % Cropping
                H.imgProcesses{end+1} = @(img)imcrop(img,getPosition(H.rectROI)) ; % Cropping
            % Colors channels
                if H.nColors>1
                    usedChannels = logical([H.redChannel.Value H.greenChannel.Value H.blueChannel.Value]) ;
                    H.imgProcesses{end+1} = @(img)img(:,:,usedChannels) ;
                    H.imgProcesses{end+1} = @(img)cat(3,img,zeros([size(img(:,:,1)) 3-size(img,3)],class(img))) ;
                end
            % To black and white ?
                if H.nColors>1 && H.toBlackAndWhite.Value
                    H.imgProcesses{end+1} = @(img)sum(img/sum(usedChannels),3,'native') ;
                end
            % Frame Normalization
                if H.normalizeEach.Value
                    H.imgProcesses{end+1} = @(img) img-min(img(:)) ;
                    H.imgProcesses{end+1} = @(img) img*(double(max(getrangefromclass(img)))/double(max(img(:)))) ;
                end
        % IMAGE DISPLAY
            % Update Image
                if H.currentFrame~=fr
                    H.currentFrame = fr ;
                    H.currentImg = H.loadFrame(H.currentFrame) ;
                    H.frameText.String = ['Current Frame: ',num2str(H.currentFrame),'/',num2str(H.nFrames)] ;
                end
            % Apply to the current image for display (without cropping)
                H.processImg = H.currentImg ;
                for p = 2:length(H.imgProcesses)
                    H.processImg = H.imgProcesses{p}(H.processImg) ;
                end
            % Display Image
                if size(H.processImg,3)==1 
                    H.ImgH.CData = repmat(H.processImg(:,:,1),[1 1 3]) ;
                else
                    H.ImgH.CData = H.processImg ;
                end
            % Finally, Crop
                H.processImg = H.imgProcesses{1}(H.processImg) ;
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



% FIGURE INITIALIZATION
    function initOptionsFigure
        % Parameters
            imgRelSize = 0.7 ; % Relative size of the image to the screen
            uiHeight = 0.05 ; % Relative vertical size of the top slider to the figure
            menuSizeLeft = 0.17 ; % Relative horizontal size of the left menu to the figure
            margin = 0.005 ;
            fontSize = 12 ;
        % Init Figure
            H.fig = figure ;
                H.fig.MenuBar = 'none' ;
                H.fig.ToolBar = 'figure' ;
                H.fig.NumberTitle = 'off' ;
                H.fig.Name = 'LOAD FRAMES: CHOOSE THE PROCESSING PARAMETERS' ;
                % Positionning
                    ratios = H.vidRes./H.fig.Position(3:4) ;
                    dimFig = imgRelSize/max(ratios)*H.vidRes ;
                    menuSizeLeft = menuSizeLeft*max(dimFig)/dimFig(1) ;
                H.fig.Position = [H.fig.Position(1:2)+H.fig.Position(3:4)/2 0 0] ...
                                    + [-1/2 0 1 0]*dimFig(1)*1/(1-menuSizeLeft) ...
                                    + [0 -1/2 0 1]*dimFig(2)*1/(1-uiHeight) ;
        % Axes containing the image preview
            H.axImg = axes('outerposition',[menuSizeLeft 0 1-menuSizeLeft 1-uiHeight]) ;
                H.axImg.LooseInset = [1 1 1 1]*0.005 ;
                H.axImg.YDir = 'reverse' ;
                H.axImg.XTick = [] ; H.axImg.YTick = [] ;
            H.ImgH = image(H.currentImg) ;
                axis(H.axImg,'tight')
                axis(H.axImg ,'equal')
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
                                            ,'FontSize',fontSize ...
                                            ... ,'HorizontalAlignment','right' ...
                                        ) ; 
        % Start|End Frames
            H.startBtn = uicontrol(H.fig,'style','pushbutton','units','normalized'...
                                            ,'position',[margin 1-2*(uiHeight+margin) menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Set Start Frame' ...
                                            ,'FontSize',fontSize ...
                                            ,'callback',@(src,evt)Global_Callback('SetStart') ...
                                        ) ; 
            H.endBtn = uicontrol(H.fig,'style','pushbutton','units','normalized'...
                                            ,'position',[margin 1-3*(uiHeight+margin) menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Set End Frame' ...
                                            ,'FontSize',fontSize ...
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
                                            ,'FontSize',fontSize ...
                                        ) ; 
        % Color Mix
            if H.nColors>1
                enableChannelsChoice = 'on' ;
            else
                enableChannelsChoice = 'off' ;
            end
            H.redChannel = uicontrol(H.fig,'style','checkbox','units','normalized'...
                                            ,'position',[margin 1-6*(uiHeight+margin) menuSizeLeft/3-2*margin uiHeight-2*margin] ...
                                            ,'string','Red' ...
                                            ,'FontSize',fontSize ...
                                            ,'Enable',enableChannelsChoice ...
                                            ,'Value',true ...
                                            ,'callback',@(src,evt)Global_Callback('normal') ...
                                        ) ; 
            H.greenChannel = uicontrol(H.fig,'style','checkbox','units','normalized'...
                                            ,'position',[margin + menuSizeLeft/3 1-6*(uiHeight+margin) menuSizeLeft/3-2*margin uiHeight-2*margin] ...
                                            ,'string','Green' ...
                                            ,'FontSize',fontSize ...
                                            ,'Enable',enableChannelsChoice ...
                                            ,'Value',true ...
                                            ,'callback',@(src,evt)Global_Callback('normal') ...
                                        ) ;  
            H.blueChannel = uicontrol(H.fig,'style','checkbox','units','normalized'...
                                            ,'position',[margin + 2*menuSizeLeft/3 1-6*(uiHeight+margin) menuSizeLeft/3-2*margin uiHeight-2*margin] ...
                                            ,'string','Blue' ...
                                            ,'FontSize',fontSize ...
                                            ,'Enable',enableChannelsChoice ...
                                            ,'Value',true ...
                                            ,'callback',@(src,evt)Global_Callback('normal') ...
                                        ) ;
            H.toBlackAndWhite = uicontrol(H.fig,'style','checkbox','units','normalized'...
                                            ,'position',[margin 1-7*(uiHeight+margin) menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Black and White' ...
                                            ,'FontSize',fontSize ...
                                            ,'Enable',enableChannelsChoice ...
                                            ,'Value',true ...
                                            ,'callback',@(src,evt)Global_Callback('normal') ...
                                        ) ; 
        % Normalization
            H.normalizeEach = uicontrol(H.fig,'style','checkbox','units','normalized'...
                                            ,'position',[margin 1-8*(uiHeight+margin) menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Norm. Each' ...
                                            ,'FontSize',fontSize ...
                                            ,'Value',false ...
                                            ,'callback',@(src,evt)Global_Callback('normal') ...
                                        ) ; 
            H.normalizeGlobal = uicontrol(H.fig,'style','checkbox','units','normalized'...
                                            ,'position',[margin 1-9*(uiHeight+margin) menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Norm. Global' ...
                                            ,'FontSize',fontSize ...
                                            ,'Value',true ...
                                            ,'callback',@(src,evt)Global_Callback('normal') ...
                                        ) ;  
        % Validation Button
            H.btnOK = uicontrol(H.fig,'style','togglebutton','units','normalized'...
                                            ,'position',[margin margin menuSizeLeft-2*margin uiHeight-2*margin] ...
                                            ,'string','Load Frames' ...
                                            ,'FontSize',fontSize ...
                                            ,'value',false ...
                                        ) ; 
        % Execute the callback once
            Global_Callback('normal') ;
    end
        
      



end