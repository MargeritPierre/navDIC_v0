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
        maximumFrameRate = 1000 ; % Hz
        % Graphical parameters
            infosTxtHeight = 16 ; % Pixels
            frameSliderWidth = .4 ; % relative to toolbar size
            frameEditWidth = .03 ; % relative to toolbar size
            frameBtnWidth = .01 ; % relative to toolbar size
        % Handles Saving Parameters
            canBeSaved = {'WorkDir'...
                            ,'Cameras'...
                            ,'DAQInputs'...
                            ,'TimeLine'...
                            ,'Images'...
                            ,'InputData'...
                            ,'Seeds'...
                            ,'Macros'...
                            ,'UserData'...
                            } ;
        
    % Process the varargin
        COMMAND = '' ;
        if ~isempty(varargin)
            if ischar(varargin{1})
                COMMAND = varargin ;
            end
        end
        
    % Is navDIC already running ?
        set(0,'ShowHiddenHandles','on') ; % Force handle visibility
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
        
    % STARTUP FILE
        navDIC_startup ;
    
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
                    hd.WorkDir.Path = pwd ;
                    hd.WorkDir.CommonName = 'img_%05i' ;
                    hd.WorkDir.ImagesExtension = '.tiff' ;
        % DEVICES
            % Cameras
                hd.Cameras = [] ;
            % DAQInputs
                hd.DAQInputs = [] ;
        % DATA
            initHandleData(false)
        % MACROS
            hd.Macros = [] ;
        % PROCESSING
            % DIC
                hd.Seeds = [] ;
            % Previews
                hd.Previews = navDICPlotPreview.empty ;
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

    % Methods for scripting
        hd.script = struct() ;
        hd.script.ui = struct() ;
        hd.script.ui.clickOn = @(btn)hd.ToolBar.MainMenu.(btn).MenuSelectedFcn(hd.ToolBar.MainMenu.(btn),[]) ;
        hd.script.ui.toggle = @(btn)set(hd.ToolBar.MainMenu.(btn),'Enable',~hd.ToolBar.MainMenu.(btn).Enable) ;
        hd.script.ui.update = @()cellfun(@(fcn)fcn(),{@updateMainMenu,@updateToolbar},'uni',false) ;
        
        
        
        
        
        
        
% ===================================================================================================================    
% INITIALIZATION FUNCTIONS
% ===================================================================================================================

  
% INITIALIZE/CLEAR DATA IN HANDLES
    function initHandleData(confirm)
        % If needed, answer the user to confirm
            if confirm && hd.nFrames>0
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
       % Init the frame controls
           hd.ToolBar.stopBtn = uicontrol('style','pushbutton'...
                                            ,'tooltip','Stop the process'...
                                            ,'units','normalized'...
                                            ,'visible','off'...
                                            ,'position',[1-frameSliderWidth-frameEditWidth-2*frameBtnWidth 0 frameBtnWidth 1] + 0.000*[1 1 -2 -2] ...
                                            ,'cdata',double(max(1:10,(1:10)')<10).*reshape([1 0 0],[1 1 3]) ...
                                            ,'callback',@(src,evt)set(src,'Visible','off') ...
                                            ) ;
           hd.ToolBar.refreshBtn = uicontrol('style','pushbutton'...
                                            ,'tooltip','Refresh'...
                                            ,'units','normalized'...
                                            ,'visible','on'...
                                            ,'position',[1-frameSliderWidth-frameEditWidth-frameBtnWidth 0 frameBtnWidth 1] + 0.000*[1 1 -2 -2] ...
                                            ,'cdata',double(((-5:5).^2+(-5:5)'.^2)<20).*reshape([0 0 1],[1 1 3]) ...
                                            ,'callback',@(src,evt)cellfun(@(fcn)fcn(),{@updateMainMenu,@updateToolbar},'uni',false) ...
                                            ) ;
           hd.ToolBar.frameSlider = uicontrol('style','slider'...
                                            ,'units','normalized'...
                                            ,'enable','off'...
                                            ,'position',[1-frameSliderWidth-frameEditWidth 0 frameSliderWidth 1]...
                                            ) ;
           hd.ToolBar.currentFrameEdit = uicontrol('style','edit'...
                                            ,'units','normalized'...
                                            ,'enable','on'...
                                            ,'position',[1-frameEditWidth 0 frameEditWidth 1]...
                                            ,'callback',@(src,evt)frameChangedFunction(src)...
                                            ) ;
           % Listener for continuous slider
                addlistener(hd.ToolBar.frameSlider, 'Value', 'PostSet',@(src,evt)frameChangedFunction(hd.ToolBar.frameSlider));
        % MAKE THE FIGURE VISIBLE
            drawnow ;
       % Add shortcuts buttons
            addButtons() ;
            %setMainMenuShortcuts() ;
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
           % Open a setup
                hd.ToolBar.MainMenu.openSetup = uimenu(hd.ToolBar.MainMenu.navDIC, ...
                                                        'Label','Open a Setup', ...
                                                        'callback',@(src,evt)openSetup) ;
           % Set working dir
                hd.ToolBar.MainMenu.setDir = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Set the Working Directory', ...
                                                        'callback',@(src,evt)setPath('menu')) ;
           % DATA IMPORT
                hd.ToolBar.MainMenu.import = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Import', ...
                                                        'Enable','on', ...
                                                        'Separator','on') ;
                    hd.ToolBar.MainMenu.importFrames = uimenu(hd.ToolBar.MainMenu.import,...
                                                            'Label','Frames', ...
                                                            'Enable','on') ;
                        hd.ToolBar.MainMenu.importImagesFolder = uimenu(hd.ToolBar.MainMenu.importFrames,...
                                                                'Label','Image Folder', ...
                                                                'Enable','on', ...
                                                                'callback',@(src,evt)import('ImageFolder')) ;
                        hd.ToolBar.MainMenu.importVideoFrames = uimenu(hd.ToolBar.MainMenu.importFrames,...
                                                                'Label','Video', ...
                                                                'Enable','on', ...
                                                                'callback',@(src,evt)import('Video')) ;
                    hd.ToolBar.MainMenu.importData = uimenu(hd.ToolBar.MainMenu.import,...
                                                            'Label','Data', ...
                                                            'Enable','on') ;
                        hd.ToolBar.MainMenu.importMatFile = uimenu(hd.ToolBar.MainMenu.importData,...
                                                                'Label','MAT File', ...
                                                                'Enable','on', ...
                                                                'callback',@(src,evt)import('MATFile')) ;
                        hd.ToolBar.MainMenu.importCSV = uimenu(hd.ToolBar.MainMenu.importData,...
                                                                'Label','CSV File', ...
                                                                'Enable','on', ...
                                                                'callback',@(src,evt)import('CSVFile')) ;
           % DATA EXPORT
                hd.ToolBar.MainMenu.export = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Export', ...
                                                        'Enable','on', ...
                                                        'Separator','off') ;
                    hd.ToolBar.MainMenu.exportImages = uimenu(hd.ToolBar.MainMenu.export,...
                                                            'Label','Images', ...
                                                            'Enable','on',...
                                                            'callback',@(src,evt)export('Images')) ;
                    hd.ToolBar.MainMenu.exportAnimation = uimenu(hd.ToolBar.MainMenu.export,...
                                                            'Label','Animation', ...
                                                            'Enable','on',...
                                                            'callback',@(src,evt)export('Animation')) ;
                    hd.ToolBar.MainMenu.exportData = uimenu(hd.ToolBar.MainMenu.export,...
                                                            'Label','Data', ...
                                                            'Enable','on',...
                                                            'callback',@(src,evt)export('Data')) ;
           % SAVE SETUP
                hd.ToolBar.MainMenu.saveSetup = uimenu(hd.ToolBar.MainMenu.navDIC,...
                                                        'Label','Save the Setup', ...
                                                        'Enable','on', ...
                                                        'Separator','on',...
                                                        'callback',@(src,evt)saveSetup()) ;
        
        % REAL-TIME ACQUISITION ----------------------------------------------------------
           hd.ToolBar.MainMenu.realTime = uimenu(hd.ToolBar.fig...
                                                ,'Label','Real-Time' ...
                                                ) ;
           % CAMERA CONNECTIONS
                hd.ToolBar.MainMenu.manageCameras = uimenu(hd.ToolBar.MainMenu.realTime,...
                                                        'Label','Manage Cameras', ...
                                                        'Separator','off', ...
                                                        'callback',@(src,evt)manageCameras) ;
           % DATA ACQUISITION INPUTS
                hd.ToolBar.MainMenu.manageInputs = uimenu(hd.ToolBar.MainMenu.realTime,...
                                                        'Label','Manage Inputs', ...
                                                        'callback',@(src,evt)manageInputs) ;
           % SAVING
                hd.ToolBar.MainMenu.saving = uimenu(hd.ToolBar.MainMenu.realTime,...
                                                        'Label','Saving', ...
                                                        'Enable','on', ...
                                                        'Checked','on', ...
                                                        'Separator','on') ;
                    hd.ToolBar.MainMenu.saveImages = uimenu(hd.ToolBar.MainMenu.saving,...
                                                            'Label','Images', ...
                                                            'Enable','on', ...
                                                            'Checked','off', ...
                                                            'callback',@(src,evt)setSaving(src)) ;
                    hd.ToolBar.MainMenu.saveInputs = uimenu(hd.ToolBar.MainMenu.saving,...
                                                            'Label','Input Data', ...
                                                            'Enable','on', ...
                                                            'Checked','off', ...
                                                            'callback',@(src,evt)setSaving(src)) ;
                    hd.ToolBar.MainMenu.saveSeeds = uimenu(hd.ToolBar.MainMenu.saving,...
                                                            'Label','Seed Data', ...
                                                            'Enable','on', ...
                                                            'Checked','off', ...
                                                            'callback',@(src,evt)setSaving(src)) ;
           % Save All Data
                hd.ToolBar.MainMenu.saveSetupData = uimenu(hd.ToolBar.MainMenu.realTime,...
                                                        'Label','Save all Acquired Data', ...
                                                        'Enable','off', ...
                                                        'callback',@(src,evt)saveSetupData()) ;
                                                    
           % Frame Rate
                hd.ToolBar.MainMenu.frameRate = uimenu(hd.ToolBar.MainMenu.realTime,...
                                                        'Label','Frame Rate', ...
                                                        ...'Enable','off', ...
                                                        'Separator','on') ;
                hd.ToolBar.MainMenu.frameRateEstimate = uimenu(hd.ToolBar.MainMenu.frameRate,...
                                                        'Label','Estimate', ...
                                                        ...'Enable','off', ...
                                                        'callback',@(src,evt)setFrameRate) ;
                hd.ToolBar.MainMenu.frameRateSet = uimenu(hd.ToolBar.MainMenu.frameRate,...
                                                        'Label','Set', ...
                                                        ...'Enable','off', ...
                                                        'Separator','on', ...
                                                        'callback',@(src,evt)setFrameRate([])) ;
                
           % Start and Stop navDIC
                hd.ToolBar.MainMenu.startStop = uimenu(hd.ToolBar.MainMenu.realTime,...
                                                        'Label','START', ...
                                                        ...'Enable','off', ...
                                                        'Separator','on', ...
                                                        'callback',@(src,evt)startContinuous) ;
                
           % Snapshot
                hd.ToolBar.MainMenu.singleShot = uimenu(hd.ToolBar.MainMenu.realTime,...
                                                        'Label','Take a Snapshot', ...
                                                        ...'Enable','off', ...
                                                        'callback',@(src,evt)singleShot) ;
                
           % RESET FRAMES
                hd.ToolBar.MainMenu.reset = uimenu(hd.ToolBar.MainMenu.realTime,...
                                                        'Label','RESET', ...
                                                        'Separator','on', ...
                                                        'callback',@(src,evt)initHandleData(true)) ;
                
        % DIC -------------------------------------------------
           hd.ToolBar.MainMenu.DIC = uimenu(hd.ToolBar.fig,'Label','DIC','Enable','off') ;
           % Manage DIC Seeds
                hd.ToolBar.MainMenu.manageDICZones = uimenu(hd.ToolBar.MainMenu.DIC,...
                                                        'Label','Manage DIC Seeds', ...
                                                        'callback',@(src,evt)manageDICZones) ;
           % Compute DIC
                hd.ToolBar.MainMenu.computeDIC = uimenu(hd.ToolBar.MainMenu.DIC,...
                                                        'Label','Compute DIC'...
                                                        ...,'Enable','off'...
                                                        ,'Separator','on'...
                                                        ) ;
               % Auto mode
                    hd.ToolBar.MainMenu.autoDIC = uimenu(hd.ToolBar.MainMenu.computeDIC,...
                                                            'Label','Auto', ...
                                                            'Checked','off', ...
                                                            'callback',@(src,evt)set(src,'Checked',~src.Checked) ...
                                                            ) ;
               % Only Non-Computed Zones
                    hd.ToolBar.MainMenu.computeSomeDIC = uimenu(hd.ToolBar.MainMenu.computeDIC,...
                                                            'Label','Non-Computed Zones Only', ...
                                                            'callback',@(src,evt)computeDIC) ;
               % All Zones
                    hd.ToolBar.MainMenu.computeAllDIC = uimenu(hd.ToolBar.MainMenu.computeDIC,...
                                                            'Label','All Zones', ...
                                                            'callback',@(src,evt)computeAllDIC) ;
                                                        
                                                        
        % MACROS -------------------------------------------------
           hd.ToolBar.MainMenu.macros = uimenu(hd.ToolBar.fig,'Label','Macros','Enable','on') ;
           % Manage Macros
                hd.ToolBar.MainMenu.manageMacros = uimenu(hd.ToolBar.MainMenu.macros,...
                                                        'Label','Manage Macros', ...
                                                        'callback',@(src,evt)manageSetupMacros) ;
           % Auto-Run Macros
                hd.ToolBar.MainMenu.autoRunMacros = uimenu(hd.ToolBar.MainMenu.macros,...
                                                        'Label','Auto Run', ...
                                                        'Separator','on', ...
                                                        'Enable','on', ...
                                                        'Checked','on', ...
                                                        'callback',@(src,evt)set(src,'checked',~get(src,'checked'))...
                                                        ) ;
           % Re-Run Macros
                hd.ToolBar.MainMenu.reRunMacros = uimenu(hd.ToolBar.MainMenu.macros,...
                                                        'Label','Re-Run', ...
                                                        'callback',@(src,evt)reRunMacros...
                                                        ) ;
                
        % VIEWS -------------------------------------------------
           hd.ToolBar.MainMenu.views = uimenu(hd.ToolBar.fig,'Label','Views');%,'Enable','off') ;
           % Preview a Camera
                hd.ToolBar.MainMenu.camPreview = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Image Preview', ...
                                                        ...'Enable','off', ...
                                                        'callback',@(src,evt)camPreview) ;
           % Preview a DIC Seed
                hd.ToolBar.MainMenu.seedPreview = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Seed Preview', ...
                                                        'callback',@(src,evt)previewSeed()) ;
           % Plot Preview
                hd.ToolBar.MainMenu.plotPreview = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Plot', ...
                                                        'callback',@(src,evt)plotPreview) ;
           % Slicing tool
                hd.ToolBar.MainMenu.sliceTool = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Images Slicing', ...
                                                        'Separator','on', ...
                                                        'callback',@(src,evt)sliceTool) ;
           % Auto Layout
                hd.ToolBar.MainMenu.autoLayout = uimenu(hd.ToolBar.MainMenu.views,...
                                                        'Label','Auto Layout', ...
                                                        ...'Enable','off', ...
                                                        'Separator','on', ...
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


% ADDS BUTTONS TOOLBAR
    function addButtons()
    end





% ===================================================================================================================    
% CONFIGURATION/IMPORT/EXPORT FUNCTIONS
% ===================================================================================================================
  
 
% SET THE WORKING DIRECTORY
    function setPath(src,varargin)
        % Open a dialog box if needed
            if strcmp(src,'menu')
                [file,path] = uiputfile('*','SELECT THE WORKING DIRECTORY, COMMON NAME AND IMAGE FORMAT','img_%05i.tiff') ;
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
%         % Select the folder of an existing setup
%             [path] = uigetdir('SELECT THE DIRECTORY OF A SAVED SETUP') ;
%             if path==0 ; return ; end
%         % Load the setup
%             [setup,hd] = loadSetup(hd,path) ;
%         % Set the WorkDir
%             setPath('open',setup.Path,setup.CommonName,setup.ImagesExtension) ;
        % Select the file to load
            [file,path] = uigetfile('.mat','Select the navDIC Setup File to Load',hd.WorkDir.Path) ;
            if path==0 ; return ; end
            filename = [path,file] ;
        % Load the file as a temporary structure
            newHD = load(filename) ;
            if isfield(newHD,'hd') ; newHD = newHD.hd ; end % old versions
        % Add/Replace the loaded fields to the current handles
            fieldsToLoad = fieldnames(newHD).' ;
            fieldsToLoad = intersect(canBeSaved,fieldsToLoad) ;
            for f = 1:length(fieldsToLoad)
                field = fieldsToLoad{f} ;
                hd.(field) = newHD.(field) ;
            end
        % If cameras are here but no Images, offer to load them
            for camID = 1:length(hd.Cameras)
                if ~isfield(newHD,'Images') ... No images where available at all
                    || length(hd.Images)<camID ... the list of was too short
                    || isempty(hd.Images{camID}) % the list contained no images for this camera
                    cam = hd.Cameras(camID) ;
                    msg = {['Camera: ',cam.Name ...
                            ,' | State: ',cam.CurrentState ...
                            ,' | Resolution: ',num2str(cam.VidObj.ROIPosition(3)),' x ' num2str(cam.VidObj.ROIPosition(4))...
                            ],'NO corresponding Images has been loaded, choose the Image Source.'} ;
                    answer = questdlg(msg,'Image Source','Image Folder','Video','Skip','Skip') ;
                    switch answer
                        case 'Image Folder'
                            [~,hd] = loadFrames(hd,'ImageFolder',camID) ;
                        case 'Video'
                            [~,hd] = loadFrames(hd,'Video',camID) ;
                    end
                end
            end
        % Compatibility with older versions
            for cam = 1:numel(hd.Images)
                if isnumeric(hd.Images{cam})
                    hd.Images{cam} = num2cell(hd.Images{cam},1:3) ;
                end
            end
        % Update the handles
            nFrames = max(cellfun(@(c)numel(c),hd.Images)) ;
            nFrames = [nFrames(:)' size(hd.TimeLine,1)] ;
            nFrames = [nFrames(:)' size(hd.InputData,1)] ;
            hd.nFrames = max(nFrames) ;
            hd.CurrentFrame = min(hd.nFrames,1) ;
        % Update the navDIC Interface
            updateMainMenu() ;
            updateToolbar() ;
        % Display infos
            disp(newline) ; disp('NEW SETUP HANDLES : ') ; display(hd) ;
    end


% IMPORT GENERAL DATA
    function import(dataType)
        % Depending on the data type...
            switch dataType
                case 'ImageFolder'
                    [valid,hd] = loadFrames(hd,dataType,length(hd.Cameras)+1) ;
                case 'Video'
                    [valid,hd] = loadFrames(hd,dataType,length(hd.Cameras)+1) ;
                case 'MATFile'
                    [valid,hd] = loadData(hd,dataType) ;
                case 'CSVFile'
                    [valid,hd] = loadData(hd,dataType) ;
            end
            if ~valid ; return ; end
        % Update handles and displays it
            clc ; disp('CURRENT navDIC HANDLES : ') ; display(hd) ;
        % Update the navDIC Interface
            updateMainMenu() ;
            updateToolbar() ;
    end


% EXPORT GENERAL DATA
    function export(dataType)
        % Depending on the data type...
            switch dataType
                case 'Images'
                    % Let the user choose a camera
                        [camID,valid] = selectCameras(hd,'single') ;
                        if ~valid ; return ; end
                    % Let the user choose an image file name
                        iLength = ceil(log10(numel(hd.Images{camID}))) ;
                        defaultName = [hd.WorkDir.Path,'\',hd.Cameras(camID).Name,'\img_%0' num2str(iLength) 'i',hd.WorkDir.ImagesExtension] ;
                        [file,path] = uiputfile(['*.' hd.WorkDir.ImagesExtension],...
                                        'SELECT THE COMMON IMAGE PATH', ... 
                                        defaultName ...
                                        ) ;
                        if path==0 ; return ; end
                        template = regexprep([path file],{'\\'},{'\\\\'}) ;
                    % Save All Images
                        wtbr = waitbar(0,'Writing Images...') ;
                        nFrames = hd.nFrames ; numel(hd.Images{camID}) ;
                        for fr = 1:nFrames 
                        % Update navDIC
                            hd.CurrentFrame = fr ;
                            updateToolbar() ;
                            hd = runMacros(hd,'onFrameChange') ;
                            hd = updateAllPreviews(hd) ;
                            drawnow ; 
                        % Write the frame
                            filename = num2str(fr,template) ;
                            imwrite(hd.Images{camID}{fr},filename) ;
                            wtbr = waitbar(fr/nFrames,wtbr,['Writing Images... (',num2str(fr),'/',num2str(hd.nFrames),')']) ;
                        end
                        delete(wtbr) ;
                case 'Animation'
                    makeAnimation()
                case 'Data'
            end
    end      


% MAKE AN ANIMATION FROM PREVIEW(S)
    function makeAnimation()
        % Update previews to clean the preview list
            hd = updateAllPreviews(hd) ;
        % Prepare the animation
            out = prepareAnimation(hd) ;
            if ~out.Valid ; warning('ANIMATION ABORTED') ; return ; end
            [~,~,ext] = fileparts(out.AnimationFile) ;
            isGIF = strcmp(lower(ext),'.gif') ;
        % Record the animation
            if ~isGIF
            % Create the Video Writer
                writerObj = VideoWriter(out.AnimationFile,'MPEG-4') ;
                writerObj.Quality = out.VideoQuality ;
                writerObj.FrameRate = out.FrameRate ;
                open(writerObj) ;
            end
            statusBox = warndlg({'The Animation is Recording...','Press OK or close to stop.'},'RECORDING...') ;
            for fr = out.FramesRecorded
                if ishandle(statusBox) % One can close it to stop the animation
                    % Update navDIC
                        hd.CurrentFrame = fr ;
                        updateToolbar() ;
                        hd = runMacros(hd,'onFrameChange') ;
                        hd = updateAllPreviews(hd) ;
                        drawnow ; 
                    % Get the individual frames
                        FR = {} ;
                        defaultColor = 255 ;
                        for ff = 1:length(out.figs)
                            FR{end+1} = getframe(out.figs(ff)) ;
                            FR{end} = FR{end}.cdata ;
                            switch out.stack
                                case 'vertical'
                                    FR{end}(:,end:max(out.sizes(:,2)),:) = defaultColor ;
                                case 'horizontal'
                                    FR{end}(end:max(out.sizes(:,1)),:,:) = defaultColor ;
                            end
                        end
                    % Create the image to record
                        switch out.stack
                            case 'vertical'
                                IMG = cat(1,FR{:}) ;
                            case 'horizontal'
                                IMG = cat(2,FR{:}) ;
                            case 'current'
                                pixPos = out.figPos./out.pixelRatio ;
                                pixPos(:,1:2) = pixPos(:,1:2) - min(pixPos(:,1:2),[],1) ;
                                pixPos = round(pixPos) ;
                                szImg = flip(max(pixPos(:,1:2)+pixPos(:,3:4),[],1) - min(pixPos(:,1:2),[],1)) ;
                                IMG = repmat(cast(defaultColor,class(FR{end})),[szImg 3]) ;
                                for ff = 1:length(out.figs)
                                    %ii = pixPos(ff,2)+(1:pixPos(ff,4)) ;
                                    ii = szImg(1)-pixPos(ff,2)-pixPos(ff,4)+(1:pixPos(ff,4)) ;
                                    jj = pixPos(ff,1)+(1:pixPos(ff,3)) ;
                                    IMG(ii,jj,:) = FR{ff} ;
                                end
                        end
                    % Add to the animation
                        if isGIF
                            [A,map] = rgb2ind(IMG,256) ;
                            if fr==out.FramesRecorded(1)
                                imwrite(A,map,out.AnimationFile,'gif','LoopCount',out.LoopCount,'DelayTime',1/out.FrameRate) ;
                            else
                                imwrite(A,map,out.AnimationFile,'gif',"WriteMode","append",'DelayTime',1/out.FrameRate) ;
                            end
                        else
                            writeVideo(writerObj,IMG) ;
                        end
                end
            end
            if ~isGIF ; close(writerObj) ; end
            if ishandle(statusBox) ; close(statusBox) ; end
    end


% SAVE THE SETUP
    function saveSetup()
        % Get the path and file 
            defName = [hd.WorkDir.Path,'/navDIC_HD.mat'] ;
            [file,path] = uiputfile('.mat','Save the navDIC Setup in a File',defName) ;
            if path==0 ; return ; end
            filename = [path,file] ;
        % Initialize the fields to be saved
            fields = canBeSaved ;
            fields = fields(ismember(fields,fieldnames(hd))) ;
        % Estimate the size of the handles and ask to save images
            if ~isempty(hd.Images)
                bytes = getfield(whos('hd'),'bytes') ;
                answer = questdlg(...
                    {'Do you want to include the Images ?',['Estimated Size: ',num2str(round(bytes/1e6)),' MB']}...
                    ,'Saving the Setup...','Yes','No','Cancel','No') ;
                if answer==0 ; return ; end
                if strcmp(answer,'Cancel') ; return ; end
                if strcmp(answer,'No') ; fields(ismember(fields,'Images')) = [] ; end
            end
        % Saving
            save(filename,'-struct','hd',fields{:},'-v7.3') ;
    end
        

        
        
% ===================================================================================================================    
% UPDATING FUNCTIONS
% ===================================================================================================================



% UPDATE INFOS TEXT
    function updateToolbar()
        % Is There any inputs ?
            nIn = 0 ; if ~isempty(hd.DAQInputs) ; nIn = length(hd.DAQInputs.Inputs) ; end
        % Infos String
            strInfos = [] ;
            strInfos = [strInfos,' ',num2str(length(hd.Cameras)),' Cameras'] ;
            strInfos = [strInfos,' | ',num2str(nIn),' DAQ.Inputs'] ;
            strInfos = [strInfos,' | ',num2str(numel(hd.Macros)),' Macros'] ;
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
        % Frame Edit
            hd.ToolBar.currentFrameEdit.String = num2str(hd.ToolBar.frameSlider.Value) ;
        % Enable or not
            if hd.nFrames>1 %&& strcmp(hd.TIMER.Running,'off') 
                hd.ToolBar.frameSlider.Enable = 'on' ;
                hd.ToolBar.currentFrameEdit.Enable = 'on' ;
            else
                hd.ToolBar.frameSlider.Enable = 'off' ;
                hd.ToolBar.currentFrameEdit.Enable = 'off' ;
            end
        % Menus
            updateMainMenu()
    end


% UPDATE THE MAIN MENU
    function updateMainMenu()
        % Is There any inputs ?
            nIn = 0 ; if ~isempty(hd.DAQInputs) ; nIn = length(hd.DAQInputs.Inputs) ; end
        % Debug blocks the mainmenu enable behavior
            if ~hd.debug
                % Data acquisition
                    if ~isempty(hd.Cameras) || nIn~=0 || ~isempty(hd.Macros)
                        hd.ToolBar.MainMenu.startStop.Enable = 'on' ;
                        hd.ToolBar.MainMenu.singleShot.Enable = 'on' ;
                        hd.ToolBar.MainMenu.frameRate.Enable = 'on' ;
                    else
                        hd.ToolBar.MainMenu.startStop.Enable = 'off' ;
                        hd.ToolBar.MainMenu.singleShot.Enable = 'off' ;
                        hd.ToolBar.MainMenu.frameRate.Enable = 'off' ;
                    end
                % Cameras
                    if ~isempty(hd.Cameras)
                        hd.ToolBar.MainMenu.camPreview.Enable = 'on' ;
                        hd.ToolBar.MainMenu.saveImages.Enable = 'on' ;
                        hd.ToolBar.MainMenu.DIC.Enable = 'on' ;
                    else
                        hd.ToolBar.MainMenu.camPreview.Enable = 'off' ;
                        hd.ToolBar.MainMenu.saveImages.Enable = 'off' ;
                        hd.ToolBar.MainMenu.DIC.Enable = 'off' ;
                    end
                % Inputs
                    if nIn~=0
                        hd.ToolBar.MainMenu.saveInputs.Enable = 'on' ;
                    else
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


% KEY PRESSED FUNCTION
    function keyPressed(src,evt)
        disp('Key Pressed !')
        disp(evt.Key)
    end

% CHANGE THE FRAME FOR PREVIEW
    function frameChangedFunction(src)
    % Format to a valid frame id
        switch src.Style
            case 'slider'
                currentFrame = src.Value ;
            case 'edit'
                currentFrame = str2double(src.String) ;
        end
        currentFrame = round(currentFrame) ;
    % Set the UIs
        hd.ToolBar.frameSlider.Value = currentFrame ;
        hd.ToolBar.currentFrameEdit.String = num2str(currentFrame) ;
    % Set navDIC state
        if hd.CurrentFrame==currentFrame ; return ; end
        hd.CurrentFrame = currentFrame ;
    % Update the navDIC state
        hd = runMacros(hd,'onFrameChange') ;
        hd = updateAllPreviews(hd) ;
        updateToolbar() ;
    end




% ===================================================================================================================    
% REAL-TIME TIMER FUNCTION
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
            % Execute macros
                if ~isempty(hd.Macros)
                    disp(['      Macros : ' num2str(t*1000,'%.1f'),' ms']) ;
                    hd = runMacros(hd,'onNewFrame') ;
                    t = toc(startTime)-lastTime ;
                    disp(['      ------ ' num2str(t*1000,'%.1f'),' ms']) ;
                    lastTime = toc(startTime) ;
                end
            % Save Acquired Data
                hd = saveCurrentFrame(hd) ;
                t = toc(startTime)-lastTime ;
                disp([' - Save : ' num2str(t*1000,'%.1f'),' ms']) ;
                lastTime = toc(startTime) ;
            % Processing
                % DIC
                if strcmp(hd.ToolBar.MainMenu.autoDIC.Checked,'on')
                    hd = updateDIC(hd) ;
                    t = toc(startTime)-lastTime ;
                    disp([' - Compute DIC : ' num2str(t*1000,'%.1f'),' ms']) ;
                    lastTime = toc(startTime) ;
                end
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


% SET THE CAMERAS
    function manageCameras
        % Stop all cameras
            hd = stopAllCameras(hd) ;
        % Open the manageMultiCameras Tool
            [hd.Cameras,camsHasChanged] = manageMultiCameras(hd.Cameras) ;
        % Re-start all cameras
            %hd = startAllCameras(hd) ;
        % If nothing changed...
            if ~camsHasChanged ; return ; end
        % Ask to clear the data
            initHandleData(true) ;
        % Update Infos
            updateToolbar() ;
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

% SET THE EXTERNAL MACROS
    function manageSetupMacros
        % Open the manageDAQInputs Tool
            [hd.Macros,macrosHasChanged] = manageMacros(hd.Macros) ;
            if ~macrosHasChanged ; return ; end
        % Ask to clear the data
            initHandleData(true) ;
        % Update Infos
            updateToolbar() ;
    end

% RE-RUN THE EXTERNAL MACROS
    function reRunMacros
        disp('Re-Run Macros')
        startFr = hd.CurrentFrame ; 1 ;
        endFr = hd.nFrames ; 1 ;
        hd.ToolBar.stopBtn.Visible = 'on' ;
        for fr = startFr:endFr
            if strcmp(hd.ToolBar.stopBtn.Visible,'off') ; break ; end
            hd.CurrentFrame = fr ;
            hd = runMacros(hd,'onFrameChange') ;
            hd = runMacros(hd,'run') ;
            hd = updateAllPreviews(hd) ;
            updateToolbar() ;
            drawnow ;
        end
        hd.ToolBar.stopBtn.Visible = 'off' ;
    end


% START CONTINUOUS SETUP
    function startContinuous()
        hd.ToolBar.frameSlider.Enable = 'off' ;
        hd.ToolBar.currentFrameEdit.Enable = 'off' ;
        hd.ToolBar.MainMenu.startStop.Label = 'STOP' ;
        hd.ToolBar.MainMenu.startStop.Callback = @(src,evt)stopContinuous ;
        hd.ToolBar.stopBtn.Visible = 'on' ;
        hd.ToolBar.stopBtn.Callback = @(src,evt)stopContinuous ;
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
        hd.ToolBar.stopBtn.Visible = 'off' ;
        hd.ToolBar.stopBtn.Callback = @(src,evt)set(src,'Visible','off') ;
    end


% TAKE A SINGLE SHOT AND PROCESS
    function singleShot()
        if strcmp(hd.TIMER.Running,'on') ; return ; end
        timerFunction() ;     
    end


% SET THE FRAME RATE
    function setFrameRate(frameRate)
    % Evaluate the frame rate ?
        if nargin==0 ; frameRate = evalMaxFrameRate() ; end
    % Ask the user to set the frame rate if needed
        if isempty(frameRate)
            frameRate = inputdlg('Set the Frame Rate (Hz)','navDIC Frame Rate',1,{num2str(hd.FrameRate)}) ;
            if isempty(frameRate) ; return ; end
            frameRate = str2num(frameRate{1}) ;
            if frameRate>maximumFrameRate 
                errordlg(['Maximum Frame Rate is ',num2str(maximumFrameRate),' Hz']) ;
                return ;
            end
        end
    % Change the timer period
        running = strcmp(hd.TIMER.Running,'on') ;
        if running ; stop(hd.TIMER) ; end
        while strcmp(hd.TIMER.Running,'on') ; end % wait for the trigger to stop
        hd.TIMER.Period = round(1/frameRate*1000)/1000 ; % millisecond precision
        if running ; start(hd.TIMER) ; end
    % Set the framerate
        hd.FrameRate = 1/hd.TIMER.Period ;
    end


% EVALUATE THE MAXIMUM FRAME RATE
    function avisedFR = evalMaxFrameRate()
        % Evaluate the maximumFrameRate by iterating the global timerFunction
            evalTime = 3 ; % seconds
            maxIt = 100 ; % maximum number of iterations
        % Backup the config
            hd_Bkp = hd ;
        % Stop the timer
            stopContinuous() ;
        % Execute it while it last less than evalTime
            itTimes = NaN(1,maxIt) ;
            it = 0 ;
            startTime = tic ;
            while it<maxIt && toc(startTime)<evalTime
                t = tic ;
                timerFunction() ;
                it = it+1 ;
                itTimes(it) = toc(t) ;
            end
        % Evaluate the maxFrameRate
            medFR = 1/median(itTimes(1:it)) ;
            avisedFR = min(0.8*medFR,maximumFrameRate) ;
        % Reset all data OK
            hd = hd_Bkp ;
        % Update toolbar and previews
            updateToolbar() ;
            hd = updateAllPreviews(hd) ;
        % Prompt the maxFrameRate
            answer = questdlg({['The Median Frame Rate is ',num2str(medFR,'%.2f'),' Hz'],...
                                ['Set the Frame Rate to ',num2str(avisedFR,'%.2f'),' Hz ?']},'Evaluated Frame rate','Yes','No','No') ;
            if strcmp(answer,'No')
                avisedFR = [] ;
            end
    end


% SET THE REAL-TIME SAVING OPTIONS
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

% SAVE THE SETUP DATA
    function saveSetupData()
        hd = saveAllFrames(hd) ;
    end




% ===================================================================================================================    
% DIGITAL IMAGE CORRELATION FUNCTIONS
% ===================================================================================================================


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
            hd.Previews(end+1) = prev ;
        end
    end

% COMPUTE NON-COMPUTED DIC ZONES
    function computeDIC
        disp('computeDIC')
    end


% COMPUTE ALL DIC ZONES
    function computeAllDIC
        disp('computeAllDIC')
        if ~any([hd.Seeds.compDisp]) ; return ; end
    % Backup DIC data
        for seed = hd.Seeds(:)'
            seed.MovingPoints_bkp = seed.MovingPoints ;
        end
    % Run DIC
        startFr = max(min([hd.Seeds([hd.Seeds.compDisp]).RefFrame]),1) ;
        endFr = min(max([hd.Seeds.FrameRange]),hd.nFrames) ;
        hd.ToolBar.stopBtn.Visible = 'on' ;
        for fr = startFr:endFr
            if strcmp(hd.ToolBar.stopBtn.Visible,'off') ; break ; end
            hd.CurrentFrame = fr ;
            hd = runMacros(hd,'onFrameChange') ;
            hd = updateDIC(hd) ;
            hd = updateAllPreviews(hd) ;
            updateToolbar() ;
            drawnow ;
        end
        hd.ToolBar.stopBtn.Visible = 'off' ;
    end



% ===================================================================================================================    
% VIEWS-RELATED FUNCTIONS
% ===================================================================================================================



% PREVIEW A CAMERA
    function camPreview
        prev = navDICCameraPreview(hd) ;
        if prev.isValid
            hd.Previews(end+1) = prev ;
        end
    end


% ADD A PLOT PREVIEW
    function plotPreview()
        hd.Previews(end+1) = navDICPlotPreview(hd) ;
    end


% START THE SLICING TOOL
    function sliceTool()
        hd.Previews(end+1) = navDICSlicingTool(hd) ;
        hd.Previews(end+1) = hd.Previews(end).SlicePrev ;
    end

% AUTO LAYOUT OF PREVIEWS
    function autoLayout 
    end



        
% ===================================================================================================================    
% navDIC INTERFACE FUNCTIONS
% ===================================================================================================================



% SWITCH IN DEBUG MODE
    function debugMode()
        hd.debug = true ;
        % Enable all Menus
            set(findobj(hd.ToolBar.fig,'type','uimenu'),'enable','on') ;
    end


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
            figsToClose = [hd.ToolBar.fig,findobj(0,'tag',navDICTag)'] ;
            set(figsToClose,'CloseRequestFcn','closereq')
            close(figsToClose) ;
    end









        
% ===================================================================================================================    
% END OF THE MAIN SCRIPT
% ===================================================================================================================

end



