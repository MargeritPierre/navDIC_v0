function navDIC(varargin)
% varargin allows to pass some special commands
% 'exit': close all figures and clear hd
% 'recover': try to recover the navDIC interface (with an existing hd)
        
% ===================================================================================================================    
% MAIN FUNCTION
% ===================================================================================================================
    
    % DEFAULT VALUES
        global hd
        navDICTag = 'navDIC' ;
        defaultFrameRate = 1 ; % Hz
        maximumFrameRate = 25 ; % Hz
        % Graphical parameters
            infosTxtHeight = 16 ; % Pixels
            frameSliderWidth = .4 ; % relative to toolbar size
        
    % Process the varargin
        COMMAND = '' ;
        if ~isempty(varargin)
            if ischar(varargin{1})
                COMMAND = varargin ;
            end
        end
        
    % Is navDIC already running ?
        set(0,'ShowHiddenHandles','on') ; % Force hadle visibility
        navDICFigs = findobj(groot,'tag',navDICTag) ;
        % If YES, prompt figures to Foreground
            if strcmp(COMMAND,'') && ~isempty(navDICFigs)
                navDICOnTop() ;
                return ;
            end
        set(0,'ShowHiddenHandles','off') ;
        
    % EXIT ?
        if strcmp(COMMAND,'exit')
            % Close all figures
                set(navDICFigs,'closerequestfcn','closereq')
                close(navDICFigs,'force')
            % Clear handles
                hd = [] ;
                clear hd
            % Exit
                return ;
        end
        
    % RECOVER ?
        if isempty(hd) && strcmp(COMMAND,'recover') % Check if there is HD to recover...
            warning('NO navDIC CONFIG TO RECOVER. Exiting...')
            return ;
        end
        if ~isempty(hd) % IF there is HD structure
            recover = true ;
            if ~strcmp(COMMAND,'recover') % The user has not asked, but...
                [order] = questdlg({'There is already a navDIC config available !','What do you want to do ?'}...
                                ,'WARNING !','Recover','Clear','Cancel','Recover') ;
                switch order
                    case 'Cancel'
                        return ;
                    case 'Clear'
                        recover = false ;
                    case 'Recover'
                        recover = true ;
                end
            end
            if recover % The user wants to recover a preceding session
                initToolBar() ;
                updateToolbar() ;
                return ;
            end
        end
    
    % Initialization of handles 
    % (type "global hd" in cmd to remote debugging access)
        hd = [] ; % Shared handles
        hd.initCompleted = false ;
        hd.navDICTag = navDICTag ;
        hd.debug = false ;
        hd.FrameRate = defaultFrameRate ;
        % Directories
            % navDic folder
                [hd.RootPath,~,~] = fileparts(which('navDIC.m')) ;
            % Working Directory
                hd.WorkDir = [] ;
        % SOURCES
            hd.Sources = [] ;
        % DEVICES
            % Cameras
                hd.Cameras = [] ;
            % DAQInputs
                hd.DAQInputs = [] ;
        % DATA
            initHandleData(false)
        % PROCESSING
            % DIC
                hd.Seeds = [] ;
            % Previews
                hd.Previews = [] ;
        % TIMER
            hd.TIMER = timer('ExecutionMode','FixedRate'...
                            ,'Period',1/defaultFrameRate ...
                            ,'TimerFcn',@(src,evt)timerFunction()...
                            ) ;
        % Initialization completed
            hd.initCompleted = true ;
        
    % No changes since last Save
        hd.hasChangedSinceLastSave = false ;
        
    % Init the Main Toolbar
        initToolBar() ;
        
        
        

        
        
% ===================================================================================================================    
% UTIL FUNCTIONS
% ===================================================================================================================

% FUNCTION EXECUTED BY THE TIMER TO TAKE A SHOT
    function timerFunction()
        % Set toolbar in red
            hd.ToolBar.infosTxt.BackgroundColor = [1 0.3 0.3] ;
            drawnow ;
        % MAIN FUNCTION
            startTime = tic ;
            disp('----- Timer Fcn -----') ;
            % Add a frame
                hd.nFrames = hd.nFrames+1 ;
                hd.CurrentFrame = hd.nFrames ;
            % Capture
                lastTime = toc(startTime) ;
                disp(' - Get Data') ;
                % Time
                    hd.TimeLine(hd.nFrames,:) = clock() ;
                    t = toc(startTime)-lastTime ;
                    disp(['      Clock : ' num2str(t*1000,'%.1f'),' ms']) ;
                    lastTime = toc(startTime) ;
                % Images
                    hd = captureCameras(hd) ;
                    t = toc(startTime)-lastTime ;
                    disp(['      Cameras : ' num2str(t*1000,'%.1f'),' ms']) ;
                    lastTime = toc(startTime) ;
                % Inputs
                    hd = captureInputs(hd) ;
                    t = toc(startTime)-lastTime ;
                    disp(['      Inputs : ' num2str(t*1000,'%.1f'),' ms']) ;
                    lastTime = toc(startTime) ;
            % Save Acquired Data
                hd = saveCurrentSetup(hd) ;
                t = toc(startTime)-lastTime ;
                disp([' - Save : ' num2str(t*1000,'%.1f'),' ms']) ;
                lastTime = toc(startTime) ;
            % Processing
                % DIC
                    hd = updateDIC(hd) ;
                    t = toc(startTime)-lastTime ;
                    disp([' - Compute DIC : ' num2str(t*1000,'%.1f'),' ms']) ;
                    lastTime = toc(startTime) ;
            % Previews
                hd = updateAllPreviews(hd) ;
                t = toc(startTime)-lastTime ;
                disp([' - Update Previews : ' num2str(t*1000,'%.1f'),' ms']) ;
                lastTime = toc(startTime) ;
        % Update Infos
            %pause(0.005) ;
            updateToolbar() ;
            t = toc(startTime)-lastTime ;
            disp([' - Update Toolbar : ' num2str(t*1000,'%.1f'),' ms']) ;
            lastTime = toc(startTime) ;
        % Time
            disp(['Total Time : ' num2str(toc(startTime)*1000,'%.1f'),' ms']) ;
            disp('------------------------') ;
            disp('') ;
        % Reset toolbar Color
            hd.ToolBar.infosTxt.BackgroundColor = hd.ToolBar.infosTxt.UserData.DefaultBackgroundColor ;
    end

% CLEAR DATA IN HANDLES
    function initHandleData(confirm)
        % If needed, answer the user to confirm
            if confirm
                answer = questdlg('DO YOU WANT TO CLEAR THE DATA ACQUIRED PREVIOUSLY ?','Clearing Data','Yes','No','No') ;
                if isempty(answer) ; return ; end
                if strcmp(answer,'No') ; return ; end
            end
        % RESET
            % Frames
                hd.nFrames = 0 ;
                hd.CurrentFrame = 0 ;
            % Time
                hd.TimeLine = [] ;
            % Images
                hd.Images = {} ;
            % Inputs
                hd.InputData = [] ;
        % Update Infos
            if hd.initCompleted 
                updateToolbar() ; 
                hd = updateAllPreviews(hd) ;
            end
    end

% UPDATE INFOS TEXT
    function updateToolbar()
        % Is There any inputs ?
            nIn = 0 ; if ~isempty(hd.DAQInputs) ; nIn = length(hd.DAQInputs.Inputs) ; end
        % Infos String
            strInfos = [] ;
            strInfos = [strInfos,' ',num2str(length(hd.Cameras)),' Cameras'] ;
            strInfos = [strInfos,' | ',num2str(nIn),' DAQ.Inputs'] ;
            strInfos = [strInfos,' | ',num2str(length(hd.Seeds)),' DIC.Seeds'] ;
            strInfos = [strInfos,' | ',num2str(length(hd.Previews)),' Previews'] ;
            strInfos = [strInfos,' | Frame ',num2str(hd.CurrentFrame),'/',num2str(hd.nFrames)] ;
            hd.ToolBar.infosTxt.String = strInfos ;
        % Frame Slider
            minVal = min(1,hd.nFrames) ;
            maxVal = max(hd.nFrames,1) ;
            hd.ToolBar.frameSlider.Min = minVal ;
            hd.ToolBar.frameSlider.Max = maxVal ;
            hd.ToolBar.frameSlider.Value = min([max(minVal,hd.CurrentFrame),maxVal,hd.nFrames]) ;
            minSliderStep = 1/max(2,hd.nFrames) ;
            maxSliderStep = max(minSliderStep,1/10) ;
            hd.ToolBar.frameSlider.SliderStep = [minSliderStep maxSliderStep] ;
            % Enable or not
                if hd.nFrames>1 %&& strcmp(hd.TIMER.Running,'off') 
                    hd.ToolBar.frameSlider.Enable = 'on' ;
                else
                    hd.ToolBar.frameSlider.Enable = 'off' ;
                end
        % Menus
            updateMainMenu()
    end

% CHANGE THE FRAME FOR PREVIEW
    function changeFrame()
        hd.ToolBar.frameSlider.Value = round(hd.ToolBar.frameSlider.Value) ;
        if hd.CurrentFrame==hd.ToolBar.frameSlider.Value ; return ; end
        hd.CurrentFrame = hd.ToolBar.frameSlider.Value ;
        hd = updateAllPreviews(hd) ;
        updateToolbar() ;
    end

% SWITCH IN DEBUG MODE
    function debugMode()
        hd.debug = true ;
        % Enable all Menus
            set(findobj(hd.ToolBar.fig,'type','uimenu'),'enable','on') ;
    end

% SET THE WORKING DIRECTORY
    function setPath(src,varargin)
        % Open a dialog box if needed
            if strcmp(src,'menu')
                [file,path] = uiputfile('*','SELECT THE WORKING DIRECTORY, COMMON NAME AND IMAGE FORMAT','img.tif') ;
                if file ==0 ; return ; end
                [~,file,ext] = fileparts(file) ;
                varargin = {path,file,ext} ;
            end
        % Update Work. Dir. Infos
            hd.WorkDir = [] ;
            hd.WorkDir.Path = varargin{1} ;
            hd.WorkDir.CommonName = varargin{2} ;
            hd.WorkDir.ImagesExtension = varargin{3} ;
        % Display infos
            disp('WORKING DIRECTORY : ') ;
            disp(hd.WorkDir) ;
    end

% OPEN A SETUP
    function openSetup
        % Select the folder of an existing setup
            [path] = uigetdir('SELECT THE DIRECTORY OF A SAVED SETUP') ;
            if path==0 ; return ; end
        % Load the setup
            [setup,hd] = loadSetup(hd,path) ;
        % Set the WorkDir
            setPath('open',setup.Path,setup.CommonName,setup.ImagesExtension) ;
        % Update handles and displays it
            clc ; disp('CURRENT SETUP HANDLES : ') ; display(hd) ;
        % Update the navDIC Interface
            updateMainMenu() ;
            updateToolbar() ;
    end

% OPEN A SETUP FROM A VIDEO
    function loadFromVideo
        % Select the video File
            [file,path] = uigetfile('*','SELECT THE VIDEO FILE TO IMPORT') ;
            if path==0 ; return ; end
        % Load the setup
            [setup,hd] = loadVideoFrames(hd,[path file]) ;
            if ~setup.Valid ; return ; end
        % Set the WorkDir
            setPath('open',setup.Path,setup.CommonName,setup.ImagesExtension) ;
        % Update handles and displays it
            clc ; disp('CURRENT SETUP HANDLES : ') ; display(hd) ;
        % Update the navDIC Interface
            updateMainMenu() ;
            updateToolbar() ;
    end

% SAVE THE SETUP DATA
    function saveSetupData()
        hd = saveAllSetupData(hd) ;
    end

% MAKE AN ANIMATION FROM PREVIEW(S)
    function makeAnimation()
        % Update previews to clean the preview list
            hd = updateAllPreviews(hd) ;
        % Prepare the animation
            out = prepareAnimation(hd) ;
            if ~out.Valid ; warning('ANIMATION ABORTED') ; return ; end
        % Record the animation
            open(out.writerObj) ;
            statusBox = warndlg({'The Animation is Recording...','Press OK or close to stop.'},'RECORDING...') ;
            for fr = out.params.FramesRecorded
                if ishandle(statusBox) % One can close it to stop the animation
                    hd.CurrentFrame = fr ;
                    updateToolbar() ;
                    hd = updateAllPreviews(hd) ;
                    drawnow ; 
                    writeVideo(out.writerObj,getframe(out.fig)) ;
                end
            end
            close(out.writerObj) ;
            if ishandle(statusBox) ; close(statusBox) ; end
    end


% START CONTINUOUS SETUP
    function startContinuous()
        hd.ToolBar.frameSlider.Enable = 'off' ;
        hd.ToolBar.MainMenu.startStop.Label = 'STOP' ;
        hd.ToolBar.MainMenu.startStop.Callback = @(src,evt)stopContinuous ;
        start(hd.TIMER) ;
    end

% STOP CONTINUOUS SETUP
    function stopContinuous() 
        if strcmp(hd.TIMER.Running,'on')
            stop(hd.TIMER) ;
            while strcmp(hd.TIMER.Running,'on') ; end
        end  
        hd.ToolBar.MainMenu.startStop.Label = 'START' ;
        hd.ToolBar.MainMenu.startStop.Callback = @(src,evt)startContinuous ;
    end

% TAKE A SINGLE SHOT AND PROCESS
    function singleShot()
        if strcmp(hd.TIMER.Running,'on') ; return ; end
        timerFunction() ;     
    end


% SET THE CAMERAS
    function manageCameras
        % Stop all cameras
            hd = stopAllCameras(hd) ;
        % Open the manageMultiCameras Tool
            [hd.Cameras,camsHasChanged] = manageMultiCameras(hd.Cameras) ;
        % Re-start all cameras
            hd = startAllCameras(hd) ;
        % If nothing changed...
            if ~camsHasChanged ; return ; end
        % Ask to clear the data
            initHandleData(true) ;
        % Update Infos
            updateToolbar() ;
    end

% PREVIEW A CAMERA
    function camPreview
        prev = navDICCameraPreview(hd) ;
        if prev.isValid
            hd.Previews{end+1} = prev ;
        end
    end

% SET THE EXTERNAL INPUTS
    function manageInputs
        % Open the manageDAQInputs Tool
            [hd.DAQInputs,inputsHasChanged] = manageDAQInputs(hd.DAQInputs) ;
            if ~inputsHasChanged ; return ; end
        % Ask to clear the data
            initHandleData(true) ;
        % Update Infos
            updateToolbar() ;
    end

% SET THE EXTERNAL INPUTS
    function inputPreview
    end

% SET THE FRAME RATE
    function setFrameRate()
        frameRate = evalMaxFrameRate() ;
        if isempty(frameRate)
            frameRate = inputdlg('Set the Frame Rate (Hz)','navDIC Frame Rate',1,{num2str(hd.FrameRate)}) ;
            if isempty(frameRate) ; return ; end
            frameRate = str2num(frameRate{1}) ;
            if frameRate>maximumFrameRate 
                errordlg(['Maximum Frame Rate is ',num2str(maximumFrameRate),' Hz']) ;
                return ;
            end
        end
        hd.FrameRate = frameRate ;
        hd.TIMER.Period = round(1/frameRate*1000)/1000 ; % millisecond precision
    end

% EVALUATE THE MAXIMUM FRAME RATE
    function avisedFR = evalMaxFrameRate()
        % Evaluate the maximumFrameRate by iterating the global timerFunction
            evalTime = 3 ; % seconds
        % Backup the config
            hd_Bkp = hd ;
        % Stop the timer
            stopContinuous() ;
        % Execute it while it last less than evalTime
            startTime = tic ;
            maxItTime = 0 ;
            it = 0 ;
            while toc(startTime)<evalTime
                t = tic ;
                timerFunction() ;
                it = it+1 ;
                maxItTime = max(maxItTime,toc(t)) ;
            end
        % Evaluate the maxFrameRate
            maxFR = 1/maxItTime ; % maxFR = it/toc(startTime) ;
            avisedFR = min(0.8*maxFR,maximumFrameRate) ;
        % Reset all data OK
            hd = hd_Bkp ;
        % Update toolbar and previews ;
            updateToolbar() ;
            hd = updateAllPreviews(hd) ;
        % Prompt the maxFrameRate
            answer = questdlg({['The Maximum Frame Rate is ',num2str(maxFR,'%.2f'),' Hz'],...
                                ['Set the Frame Rate to ',num2str(avisedFR,'%.2f'),' Hz ?']},'Evaluated Frame rate','Yes','No','No') ;
            if strcmp(answer,'No')
                avisedFR = [] ;
            end
            
    end

% SET THE FRAME RATE
    function setSaving(src)
        switch src.Checked
            case 'off' % The user wants to activate
                switch src.Label
                    case 'Images'
                        if isempty(hd.Cameras) ; warning('No image will be saved until a camera is added') ; end
                    case 'Input Data'
                        if isempty(hd.DAQInputs) || isempty(hd.DAQInputs.Inputs) ; warning('No data will be saved until a DAQ input is added') ; end
                    case 'Seed Data'
                        if isempty(hd.Seeds) ; warning('No data will be saved until a seed is added') ; end
                end
                src.Checked = 'on' ;
            case 'on' % The user wants to desactivate
                src.Checked = 'off' ;
        end
    end

% MANAGE DIC ZONES
    function manageDICZones
        hd = manageDICSeeds(hd) ;
        % Update Infos
            updateToolbar() ;
    end

% PREVIEW AN INDIVIDUAL DIC SEED
    function previewSeed()
        if isempty(hd.Seeds) ; return ; end
        prev = navDIC2DSeedPreview(hd) ;
        if prev.isValid
            hd.Previews{end+1} = prev ;
        end
    end

% COMPUTE NON-COMPUTED DIC ZONES
    function computeDIC
        disp('computeDIC')
    end

% COMPUTE ALL DIC ZONES
    function computeAllDIC
        disp('computeAllDIC')
        for fr = 1:hd.nFrames
            hd.CurrentFrame = fr ;
            hd = updateDIC(hd) ;
            hd = updateAllPreviews(hd) ;
            updateToolbar() ;
            drawnow ;
        end
    end

% ADD A PLOT PREVIEW
    function plotPreview()
        hd.Previews{end+1} = navDICPlotPreview(hd) ;
    end

% ADD A PLOT PREVIEW
    function sliceTool()
        hd.Previews{end+1} = navDICSlicingTool(hd) ;
    end

% MANAGE AXES
    function manageViews
    end

% AUTO LAYOUT
    function autoLayout
    end



        
% ===================================================================================================================    
% GRAPHICAL FUNCTIONS
% ===================================================================================================================

% BRING ALL NAVDIC FIGURES TO FRONT
    function navDICOnTop()
        disp('tofront') ;
        navDICFigs = findobj(groot,'tag',navDICTag) ;
        if ~isempty(navDICFigs)
            for f = 1:length(navDICFigs)
                figure(navDICFigs(f)) ;
            end
        end
    end

% CREATES THE MAIN FIGURE FOR MENU AND TOOLBAR
    function initToolBar()
       % Figure creation
           screenPos = get(groot,'monitorpositions') ;
           hd.ToolBar.fig = figure('Name','navDIC v0.0',...
                                    'toolbar','none',...
                                    'menubar','none',...
                                    'outerposition',screenPos(end,:),...
                                    'dockcontrols','off',...
                                    'NumberTitle','off',...
                                    'Visible','off',...
                                    'tag',navDICTag...
                                    ) ;
           hd.ToolBar.fig.ButtonDownFcn = @(src,evt)navDICOnTop() ;
       % Close Callback
           hd.ToolBar.fig.CloseRequestFcn = @(src,evt)closeAll() ;
       % Add the main menu
           addMainMenu() ;
           updateMainMenu() ;
       % Put The toolbar at the top of the screen
           drawnow ; % I don't know why, but i'm forced to draw here...
           hd.ToolBar.fig.Position(4) = infosTxtHeight ;
           hd.ToolBar.fig.OuterPosition(2) = screenPos(end,4) + screenPos(end,2) - hd.ToolBar.fig.OuterPosition(4) ;
           drawnow ;
       % Init the Info textbox
           hd.ToolBar.infosTxt = uicontrol('style','text'...
                                            ,'fontname','Consolas'...
                                            ,'string',''...
                                            ,'units','normalized'...
                                            ,'position',[0 0 1 1]...
                                            ,'fontunits','pixels'...
                                            ,'fontsize',0.8*infosTxtHeight...
                                            ,'horizontalalignment','left'...
                                            ) ;
           hd.ToolBar.infosTxt.UserData.DefaultBackgroundColor = hd.ToolBar.infosTxt.BackgroundColor ;
       % Init the step Slider
           hd.ToolBar.frameSlider = uicontrol('style','slider'...
                                            ,'units','normalized'...
                                            ,'enable','off'...
                                            ,'position',[1-frameSliderWidth 0 frameSliderWidth 1]...
                                            ) ;
           % Listener for continuous slider
                addlistener(hd.ToolBar.frameSlider, 'Value', 'PostSet',@(src,evt)changeFrame());
        % MAKE THE FIGURE VISIBLE
            drawnow ;
       % Add shortcuts buttons
            addButtons() ;
            setMainMenuShortcuts() ;
            hd.ToolBar.fig.KeyPressFcn = @(src,evt)keyPressed(src,evt) ;
            drawnow ;
        % Set the figure handle invisible
            hd.ToolBar.fig.Visible = 'on' ;
            hd.ToolBar.fig.HandleVisibility = 'off' ;
    end


% BUILD THE MAIN MENU
    function addMainMenu()
        
        % NAVDIC ----------------------------------------------------------
           hd.ToolBar.MainMenu.navDIC = uimenu(hd.ToolBar.fig...
                                                ,'Label',navDICTag ...
                                                ) ;
           % Set working dir
                hd.ToolBar.MainMenu.setDir = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Set the Working Directory', ...
                                                        'callback',@(src,evt)setPath('menu')) ;
           % Open a setup
                hd.ToolBar.MainMenu.openSetup = uimenu(hd.ToolBar.MainMenu.navDIC, ...
                                                        'Label','Open a Setup', ...
                                                        'callback',@(src,evt)openSetup) ;
           % Load From a Video
                hd.ToolBar.MainMenu.loadFromVideo = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Load From Video', ...
                                                        'callback',@(src,evt)loadFromVideo, ...
                                                        'Separator','off') ;
           % SAVING
                hd.ToolBar.MainMenu.saving = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Saving', ...
                                                        'Enable','off', ...
                                                        'Checked','on', ...
                                                        'Separator','on') ;
                    hd.ToolBar.MainMenu.saveImages = uimenu(hd.ToolBar.MainMenu.saving,...
                                                            'Label','Images', ...
                                                            'Enable','off', ...
                                                            'Checked','on', ...
                                                            'callback',@(src,evt)setSaving(src)) ;
                    hd.ToolBar.MainMenu.saveInputs = uimenu(hd.ToolBar.MainMenu.saving,...
                                                            'Label','Input Data', ...
                                                            'Enable','off', ...
                                                            'Checked','on', ...
                                                            'callback',@(src,evt)setSaving(src)) ;
                    hd.ToolBar.MainMenu.saveSeeds = uimenu(hd.ToolBar.MainMenu.saving,...
                                                            'Label','Seed Data', ...
                                                            'Enable','off', ...
                                                            'Checked','on', ...
                                                            'callback',@(src,evt)setSaving(src)) ;
           % Save All Data
                hd.ToolBar.MainMenu.saveSetupData = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Save all Acquired Data', ...
                                                        'Enable','off', ...
                                                        'callback',@(src,evt)saveSetupData()) ;
           % Save All Data
                hd.ToolBar.MainMenu.exportAnimation = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Export Animation', ...
                                                        'Enable','on', ...
                                                        'callback',@(src,evt)makeAnimation()) ;
                                                    
           % Frame Rate
                hd.ToolBar.MainMenu.frameRate = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Frame Rate', ...
                                                        ...'Enable','off', ...
                                                        'Separator','on', ...
                                                        'callback',@(src,evt)setFrameRate) ;
                
           % Start and Stop navDIC
                hd.ToolBar.MainMenu.startStop = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','START', ...
                                                        ...'Enable','off', ...
                                                        'Separator','on', ...
                                                        'callback',@(src,evt)startContinuous) ;
                
           % Snapshot
                hd.ToolBar.MainMenu.singleShot = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Take a Snapshot', ...
                                                        ...'Enable','off', ...
                                                        'callback',@(src,evt)singleShot) ;
                
           % RESET FRAMES
                hd.ToolBar.MainMenu.reset = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','RESET', ...
                                                        'Separator','on', ...
                                                        'callback',@(src,evt)initHandleData(true)) ;
                
        % CAMERAS ---------------------------------------------------------
           hd.ToolBar.MainMenu.cameras = uimenu(hd.ToolBar.fig,'Label','Cameras') ;
           % Set Cameras
                hd.ToolBar.MainMenu.manageCameras = uimenu(hd.ToolBar.MainMenu.cameras,...
                                                        'Label','Manage Cameras', ...
                                                        'callback',@(src,evt)manageCameras) ;
           % Preview a Camera
                hd.ToolBar.MainMenu.camPreview = uimenu(hd.ToolBar.MainMenu.cameras,...
                                                        'Label','Preview a Camera', ...
                                                        ...'Enable','off', ...
                                                        'callback',@(src,evt)camPreview) ;
                
        % EXTERNAL INPUTS -------------------------------------------------
           hd.ToolBar.MainMenu.extInputs = uimenu(hd.ToolBar.fig,'Label','Inputs') ;
           % Set Inputs
                hd.ToolBar.MainMenu.manageInputs = uimenu(hd.ToolBar.MainMenu.extInputs,...
                                                        'Label','Manage Inputs', ...
                                                        'callback',@(src,evt)manageInputs) ;
           % Preview an Input
                hd.ToolBar.MainMenu.inputPreview = uimenu(hd.ToolBar.MainMenu.extInputs,...
                                                        'Label','Preview an Input', ...
                                                        ...'Enable','off', ...
                                                        'callback',@(src,evt)inputPreview) ;
                
        % DIC -------------------------------------------------
           hd.ToolBar.MainMenu.DIC = uimenu(hd.ToolBar.fig,'Label','DIC','Enable','off') ;
           % Manage DIC Seeds
                hd.ToolBar.MainMenu.manageDICZones = uimenu(hd.ToolBar.MainMenu.DIC,...
                                                        'Label','Manage DIC Seeds', ...
                                                        'callback',@(src,evt)manageDICZones) ;
           % Preview a DIC Seed
                hd.ToolBar.MainMenu.previewSeed = uimenu(hd.ToolBar.MainMenu.DIC,...
                                                        'Label','Preview a Seed', ...
                                                        'callback',@(src,evt)previewSeed()) ;
           % Compute DIC
                hd.ToolBar.MainMenu.computeDIC = uimenu(hd.ToolBar.MainMenu.DIC,...
                                                        'Label','Compute DIC'...
                                                        ...,'Enable','off'...
                                                        ,'Separator','on'...
                                                        ) ;
               % Only Non-Computed Zones
                    hd.ToolBar.MainMenu.computeSomeDIC = uimenu(hd.ToolBar.MainMenu.computeDIC,...
                                                            'Label','Non-Computed Zones Only', ...
                                                            'callback',@(src,evt)computeDIC) ;
               % All Zones
                    hd.ToolBar.MainMenu.computeAllDIC = uimenu(hd.ToolBar.MainMenu.computeDIC,...
                                                            'Label','All Zones', ...
                                                            'callback',@(src,evt)computeAllDIC) ;
                
        % VIEWS -------------------------------------------------
           hd.ToolBar.MainMenu.views = uimenu(hd.ToolBar.fig,'Label','Views');%,'Enable','off') ;
           % Plot Preview
                hd.ToolBar.MainMenu.plotPreview = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Plot Preview', ...
                                                        'callback',@(src,evt)plotPreview) ;
           % Slicing tool
                hd.ToolBar.MainMenu.sliceTool = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Slicing', ...
                                                        'callback',@(src,evt)sliceTool) ;
           % Manage Axes
                hd.ToolBar.MainMenu.manageViews = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Manage Views', ...
                                                        'callback',@(src,evt)manageViews) ;
           % Auto Layout
                hd.ToolBar.MainMenu.autoLayout = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Auto Layout', ...
                                                        ...'Enable','off', ...
                                                        'callback',@(src,evt)autoLayout) ;
                
        % HELP -------------------------------------------------
           hd.ToolBar.MainMenu.help = uimenu(hd.ToolBar.fig,'Label','?');%,'Enable','off') ;
           % Manage Axes
                hd.ToolBar.MainMenu.about = uimenu(hd.ToolBar.MainMenu.help,...
                                                        'Label','About', ...
                                                        'callback',@(src,evt){}) ;
           % DebugMode
                hd.ToolBar.MainMenu.debugMode = uimenu(hd.ToolBar.MainMenu.help,...
                                                        'Label','Debug', ...
                                                        ...'Enable','on', ...
                                                        'callback',@(src,evt)debugMode) ;
                                                    
                                                    
    end

% KEY PRESSED FUNCTION
    function keyPressed(src,evt)
        disp('Key Pressed !')
        disp(evt.Key)
    end

% SET THE MAIN TOOLBAR SHORTCUTS
    function setMainMenuShortcuts()
        % Define Shortcuts
            shortcuts = {'Open a Setup','shift O'; ...
                         'Set the Working Directory','shift S'; ...
                         'START','shift ENTER'; ...
                         'Take a Snapshot','shift SPACE'; ...
                         'Load From Video','shift V'; ...
                         'Frame Rate','shift F'; ...
                         'Preview a Camera','shift C'; ...
                         'Preview an Input','shift I'; ...
                         'Preview a Seed','shift D'; ...
                         'Plot Preview','shift P'; ...
                         } ;
        % Get Menus Handles
            hMenus = findobj(hd.ToolBar.fig,'type','uimenu') ; 
            menuLabels = {hMenus.Label} ;
        % Set "accelerators" (or shortcuts)
            for s = 1:size(shortcuts,1)
                idHandle = ismember(menuLabels,{shortcuts{s,1}}) ;
                if any(idHandle)
                    for it = 1:5 % Try multiple times to find the java object
                    jHandle = findjobj(hMenus(idHandle)) ;
                        if ~isempty(jHandle)
                            jAccelerator = javax.swing.KeyStroke.getKeyStroke(shortcuts{s,2}) ;
                            jHandle.setAccelerator(jAccelerator) ;
                            continue
                        end
                    end
               end
            end
    end


% UPDATE THE MAIN MENU
    function updateMainMenu()
        % Is There any inputs ?
            nIn = 0 ; if ~isempty(hd.DAQInputs) ; nIn = length(hd.DAQInputs.Inputs) ; end
        % Debug blocks the mainmenu enable behavior
            if ~hd.debug
                % Data acquisition
                    if ~isempty(hd.Cameras) || nIn~=0
                        hd.ToolBar.MainMenu.startStop.Enable = 'on' ;
                        hd.ToolBar.MainMenu.singleShot.Enable = 'on' ;
                        hd.ToolBar.MainMenu.frameRate.Enable = 'on' ;
                        hd.ToolBar.MainMenu.saveImages.Enable = 'on' ;
                    else
                        hd.ToolBar.MainMenu.startStop.Enable = 'off' ;
                        hd.ToolBar.MainMenu.singleShot.Enable = 'off' ;
                        hd.ToolBar.MainMenu.frameRate.Enable = 'off' ;
                        hd.ToolBar.MainMenu.saveImages.Enable = 'off' ;
                    end
                % Cameras
                    if ~isempty(hd.Cameras)
                        hd.ToolBar.MainMenu.camPreview.Enable = 'on' ;
                        hd.ToolBar.MainMenu.DIC.Enable = 'on' ;
                    else
                        hd.ToolBar.MainMenu.camPreview.Enable = 'off' ;
                        hd.ToolBar.MainMenu.DIC.Enable = 'off' ;
                    end
                % Inputs
                    if nIn~=0
                        hd.ToolBar.MainMenu.inputPreview.Enable = 'on' ;
                        hd.ToolBar.MainMenu.saveInputs.Enable = 'on' ;
                    else
                        hd.ToolBar.MainMenu.inputPreview.Enable = 'off' ;
                        hd.ToolBar.MainMenu.saveInputs.Enable = 'off' ;
                    end
                % Seeds
                    if ~isempty(hd.Seeds)
                        hd.ToolBar.MainMenu.seedPreview.Enable = 'on' ;
                        hd.ToolBar.MainMenu.saveSeeds.Enable = 'on' ;
                    else
                        hd.ToolBar.MainMenu.seedPreview.Enable = 'off' ;
                        hd.ToolBar.MainMenu.saveSeeds.Enable = 'off' ;
                    end
            end
    end


% ADDS BUTTONS TOOLBAR
    function addButtons()
    end


% CLOSE THE PROGRAM
    function closeAll()
        % Ask the user about LOOSING INFO
            if exist('hd','var') && hd.hasChangedSinceLastSave
                button = questdlg('Close the navDIC ? All data will be lost.','CLOSE REQUEST','Yes','No','No') ;
                switch button
                    case 'No' % Stop !
                        return ;
                    case 'Yes' % Continue
                end
            end
        % STOP THE TIMER
            stop(hd.TIMER) ;
            while strcmp(hd.TIMER.Running,'on') ; end
        % STOP ALL CAMERAS
            hd = stopAllCameras(hd) ;
        % Close all figures belonging to navDIC 
            figsToClose = [hd.ToolBar.fig,findobj(0,'tag',navDICTag)] ;
            set(figsToClose,'CloseRequestFcn','closereq')
            close(figsToClose) ;
    end

end



