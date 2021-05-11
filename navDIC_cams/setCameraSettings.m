function cam = setCameraSettings(cam)
%%% [cam,valid] = setCameraSettings(cam)
% 

    % Figure tag (used to prompt figure in first plane)
        figTag = ['setCameraSettings(',cam.Name,')'] ;
        % Is the camera already in setting mode ?
            prevFig = findobj(groot,'tag',figTag) ;
            if ~isempty(prevFig)
                figure(prevFig) ;
                return ;
            end
       
        
        
    % Init INFOS Structure
        INFOS = [] ;
        initINFOS() ;
        if ~INFOS.isValid
            disp('Camera setting cancelled') ;
            return ;
        end
        INFOS = rmfield(INFOS,'isValid') ;
        
    % FULL camera ROI
        cam.VidObj.ROIPosition = [0 0 cam.VidObj.VideoResolution] ;
        
    % Init the preview timer
        cam.CurrentState = 'setting' ;
        PREVIEW = [] ;
        
    % Init the Figure
        fig = [] ;
        % FIGURE Parameters
            figRelSize = .6 ; % Relative size to screen
            menuSize = .3 ; % Relative to image width
            histoHeight = .15 ; %relative to figure height
            marginMenu = .01 ;
            marginAxes = .01 ;
            marginUIs = .005 ;
            uiHeight = .5 ;
            nBinsHistogram = 50 ;
        % Go
            initFigure() ;

    % Wait for the user to close the figure 
        while strcmp(cam.CurrentState,'setting')
            try 
                drawnow ; 
            end
        end
        
    % Set the correct ROI
        cam.VidObj.ROIPosition = INFOS.Custom.ROI.Value ;
    
    % Return the tuned cam
        return ;
        
        
        
% ===================================================================================================================    
% CHILD FUNCTIONS
% ===================================================================================================================


% UPDATE CAMERA PREVIEW
    function updatePreview(obj,event,hImage)
        % Display the current image frame. 
            frame0 = event.Data ;
            %frame0 = double(frame0) ;
            %frame0 = frame0*(1./max(getrangefromclass(event.Data))) ;
            %frame0 = double(getsnapshot(obj)) ;
        % Processing on the frame
            switch PREVIEW.derivBtn.String{PREVIEW.derivBtn.Value}
                case 'gradient'
                    frame = abs(diff(frame0([1:end,end],:),1,1)+1i*diff(frame0(:,[1:end,end]),1,2)) ;
                case 'laplacian'
                    frame = abs(diff(frame0([1,1:end,end],:),2,1)+1i*diff(frame0(:,[1,1:end,end]),2,2)) ;
                case 'noise level'
                    frame = abs(frame0-PREVIEW.LastFrame) ;
                otherwise
                    frame = frame0 ;
            end
        % Show the Acquired Frame
            hImage.CData = frame ;
        % If a frame is available, process
            % Histogram on the ROI only
                frame = imcrop(frame,INFOS.Custom.ROI.Value) ;
                updateHisto(frame) ;
        % Update Last Frame
            PREVIEW.LastFrame = frame0 ;
    end

    function updateHisto(frame)
        % Try to fing an Histogram
            histObj = findobj(fig,'type','Histogram') ;
        % Update (or create) It
            nBands = size(frame,3) ;
            if nBands==1 % Grayscale
                if isempty(histObj)
                    histogram(PREVIEW.AxHistObj,frame(:,:,1),nBinsHistogram,'facecolor','k') ;
                else
                    histObj.Data = frame(:,:,1) ;
                end
            else % Colorscale
                colors = ['r';'g';'b'] ;
                for b = 1:nBands
                    if isempty(histObj)
                        histoPreview(b) = histogram(frame(:,:,b),nBinsHistogram,'facecolor',colors(b),'facealpha',.3) ;
                    else
                        histObj(b).Data = frame(:,:,b) ;
                    end
                end
            end
        % Axes formatting
            %maxFrame = max(frame(:)) ;
            switch class(frame)
                case 'double'
                    PREVIEW.AxHistObj.XLim = [0 1] ;
                otherwise
                    PREVIEW.AxHistObj.XLim = [0 1].*max(getrangefromclass(frame)) ;
            end
    end


% UPDATE SETTINGS
    function updateSettings(src,evt)
        % Get the source identity
            uiTag = get(src,'Tag') ;
            UI = findobj(fig,'Tag',uiTag) ;
        % Is the source a IMAQ parameter ?
            if isfield(INFOS.IMAQ,uiTag)
                % Retrieve value
                    switch UI.Style
                        case 'slider'
                            value = UI.Value ;
                        case 'edit'
                            uiInfo = INFOS.IMAQ(uiTag) ;
                            switch uiInfo.Type
                                case 'string'
                                    value = UI.String ;
                                case 'double'
                                    value = str2double(UI.String) ;
                                case 'integer'
                                    value = round(str2double(UI.String)) ;
                            end
                        case 'popupmenu'
                            value = UI.String{UI.Value} ;
                    end
                % Apply
                    %setfield(cam.VidObj.Source,uiTag,value) ;
                    set(cam.VidObj.Source,uiTag,value) ;
                % Terminate
                    return ;
            end
        % ELSE It is a custom parameter
            switch uiTag
                case 'rectROI'
                    pos = round(evt) ;
                    uiROI = findobj(fig,'tag','ROI') ;
                    uiROI.String = mat2str(pos) ;
                    INFOS.Custom.ROI.Value = pos ;
                case 'Name'
                    cam.Name = UI.String ;
            end
    end


% CREATE THE UI MENU
    function createMenu()
        % Get ALL parameters INFOS
            paramsInfos = {} ;
            paramsFieldNames = {} ;
            for infoType = fieldnames(INFOS)'
                params = getfield(INFOS,infoType{:}) ;
                paramsInfos = cat(1,struct2cell(params),paramsInfos) ;
                paramsFieldNames = cat(1,fieldnames(params),paramsFieldNames) ;
            end
            paramsInfos = cell2struct(paramsInfos,paramsFieldNames) ;
        % Create the uicontrols
            nUIControls = length(paramsFieldNames) ;
            uiTotalHeight = (1-3*marginMenu-histoHeight-(nUIControls-1)*marginUIs)/nUIControls ;
            for ui = 1:nUIControls
                uiTag = paramsFieldNames{end-ui+1} ;
                uiInfos = getfield(paramsInfos,uiTag) ;
                uiPos = [marginMenu 2*marginMenu+histoHeight+(ui-1)*(marginUIs+uiTotalHeight) menuSize/(1+menuSize)-2*marginMenu uiTotalHeight] ;
                % Main UI
                    % Init
                        theUI = uicontrol('tag',uiTag,...
                                        'Units','normalized',...
                                        'Position',uiPos.*[1 1 1 uiHeight]+[0 0 0 0],...
                                        'BackgroundColor',[1 1 1]*.95) ;
                    % Individual infos
                        switch uiInfos.Type
                            case 'integer'
                                switch uiInfos.Constraint
                                    case 'bounded'
                                        theUI.Style = 'slider' ;
                                        theUI.Min = uiInfos.ConstraintValue(1) ;
                                        theUI.Max = uiInfos.ConstraintValue(2) ;
                                    case 'none'
                                        theUI.Style = 'edit' ;
                                end
                            case 'double'
                                switch uiInfos.Constraint
                                    case 'bounded'
                                        theUI.Style = 'slider' ;
                                        theUI.Min = uiInfos.ConstraintValue(1) ;
                                        theUI.Max = uiInfos.ConstraintValue(2) ;
                                    case 'none'
                                        theUI.Style = 'edit' ;
                                end
                            case 'string'
                                switch uiInfos.Constraint
                                    case 'enum'
                                        theUI.String = uiInfos.ConstraintValue ;
                                        theUI.Style = 'popupmenu' ;
                                    case 'none'
                                        theUI.Style = 'edit' ;
                                end
                            case 'ROI'
                                theUI.Style = 'text' ;
                            otherwise
                                theUI.Style = 'edit' ;
                        end
                    % Callback
                        theUI.Callback = @(src,evt)updateSettings(src,evt) ;
                % Tag UI
                    uicontrol('style','text',...
                                'string',uiTag,...
                                'Tag','legend',...
                                'Units','normalized',...
                                'Position',uiPos.*[1 1 1 1-uiHeight]+[0 uiTotalHeight*uiHeight 0 0],...
                                'BackgroundColor',[1 1 1]*.8) ;

            end
    end


% INIT UI VALUES
    function initUIValues()
        % Get all UIs
            UIs = findobj(fig,'type','uicontrol') ;
        % Treat each UI individually
            for ui = 1:length(UIs)
                % IMAQ-related UIs
                    if isfield(INFOS.IMAQ,UIs(ui).Tag)
                        switch UIs(ui).Style
                            case 'slider'
                                %UIs(ui).Value = getfield(cam.VidObj.Source,UIs(ui).Tag) ;
                                UIs(ui).Value = getfield(cam.VidObj.Source,UIs(ui).Tag) ;
                            case 'edit'
                                %UIs(ui).String = string(getfield(cam.VidObj.Source,UIs(ui).Tag)) ;
                                UIs(ui).String = num2str(getfield(cam.VidObj.Source,UIs(ui).Tag)) ;
                            case 'popupmenu'
                                %str = string(getfield(cam.VidObj.Source,UIs(ui).Tag)) ;
                                str = num2str(getfield(cam.VidObj.Source,UIs(ui).Tag)) ;
                                [~,id] = intersect(UIs(ui).String,str) ;
                                UIs(ui).Value = id ; 
                        end
                        continue ;
                    end
                % Custom-related UIs
                    switch UIs(ui).Tag
                        case 'Name'
                            UIs(ui).String = INFOS.Custom.Name.Value ; 
                        case 'ROI'
                            UIs(ui).String = mat2str(INFOS.Custom.ROI.Value) ; 
                    end
            end 
    end
      

% INIT INFOS STRUCTURE
    function initINFOS() 
    % Get Camera Parameters
        INFOS.isValid = false ;
        % PARAMS Structure is inspired from propinfo() results structure
            infoStruct = [] ;
            infoStruct.Type = {} ;
            infoStruct.Value = {} ;
            infoStruct.Constraint = {} ;
            infoStruct.ConstraintValue = {} ;
            infoStruct.DefaultValue = {} ;
            infoStruct.ReadOnly = {} ;
            infoStruct.DeviceSpecific = {} ;
            infoStruct.Accessible = {} ;
        % IMAQ Properties
            INFOS.IMAQ = cam.Infos.IMAQ ; %propinfo(cam.VidObj.Source) ;
            % Default Parameters
                defParams = {...
                                'AutoExposure';...
                                'Brightness';...
                                'FrameTimeout';...
                                'Gain';...
                                'GainAbsolute';...
                                'GainControl';...
                                'GainMode';...
                                'Gamma';...
                                'NormalizedBytesPerPacket';...
                                'Shutter';...
                                'ShutterMode';...
                                'PacketSize';...
                                'ExposureAuto';...
                                'ExposureAutoAdjustTol';...
                                'ExposureAutoAlg';...
                                'ExposureAutoMax';...
                                'ExposureAutoMin';...
                                'ExposureAutoOutliers';...
                                'ExposureAutoRate';...
                                'ExposureAutoTarget';...
                                'ExposureMode';...
                                } ;
            % Available Parameters
                paramsFields = fieldnames(INFOS.IMAQ) ;
            % Default selected IDs
                defIDs = [] ;
                for p = 1:length(defParams)
                    defIDs = [defIDs find(ismember(paramsFields,defParams{p}))] ;
                end
            % Choose list
                [paramIDs,valid] = listdlg('PromptString','Select Cameras :',...
                                            'SelectionMode','multiple',...
                                            'initialValue',defIDs,...
                                            'ListString',paramsFields) ;
                if ~valid ; return ; end
            % Exclude some parameters
%                 excludeParams = {'Parent',...
%                                 'SourceName',...
%                                 'BytesPerPacket',...
%                                 'TriggerDelay',...
%                                 'TriggerDelayControl',...
%                                 'TriggerDelayAbsolute',...
%                                 'TriggerDelayMode',...
%                                 'TriggerParameter',...
%                                 'Selected',...
%                                 'Tag',...
%                                 'Type',...
%                                 'UniqueID'...
%                                 } ;
                excludeParams = setdiff(paramsFields,paramsFields(paramIDs)) ;
                for p = 1:length(excludeParams)
                    if isfield(INFOS.IMAQ,excludeParams{p})
                        INFOS.IMAQ = rmfield(INFOS.IMAQ,excludeParams{p}) ;
                    end
                end
        % Custom properties
            INFOS.Custom = [] ;
            % Name
                INFOS.Custom.Name = infoStruct ;
                INFOS.Custom.Name.Type = 'string' ;
                INFOS.Custom.Name.Value = cam.Name ;
                INFOS.Custom.Name.Constraint = 'none' ;
                INFOS.Custom.Name.DefaultValue = cam.Infos.DeviceName ;
            % ROI
                INFOS.Custom.ROI = infoStruct ;
                INFOS.Custom.ROI.Type = 'ROI' ;
                INFOS.Custom.ROI.Value = cam.VidObj.ROIPosition ; %cam.VidObj.ROIPosition ;
                INFOS.Custom.ROI.Constraint = 'bounded' ;
        % Valid Output
            INFOS.isValid = true ;
    end
        
            
% INIT THE FIGURE
    function initFigure()
        % Video Parameters
            vidRes = cam.VidObj.VideoResolution ;
            nBands = cam.VidObj.NumberOfBands ;
            monitPos = get(groot,'monitorpositions') ;
            monitPos = monitPos(end,:) ;
            monitSize = monitPos(3:4) ;
            monit2VidRatio = min(monitSize./vidRes) ;
        % Create Figure
            fig = figure('tag',figTag,...
                        'Name','SET CAMERA SETTINGS',...
                        ...'windowstyle','modal',...
                        'toolbar','none',...
                        'menubar','none',...
                        'NumberTitle','off'...
                        ) ;
            fig.Units = 'pixels' ;
            fig.Position = vidRes(1)*[-.5 0 1 0]*figRelSize*monit2VidRatio*(1+menuSize) ...
                            + vidRes(2)*[0 -.5 0 1]*figRelSize*monit2VidRatio ...
                            + [fig.Position(1)+fig.Position(3)/2 fig.Position(2)+fig.Position(4)/2 0 0] ;
            fig.Units = 'normalized' ;
            fig.CloseRequestFcn = @(src,evt)closeFigure ;
        % Setup the preview axes
            axPreview = axes() ;
                axPreview.Units = 'pixels' ;
                axPreview.Position = [menuSize*vidRes(1) 0 vidRes(1) vidRes(2)]*figRelSize*monit2VidRatio + [1 1 0 0] ;
                axPreview.Units = 'normalized' ;
                axPreview.Position = axPreview.Position + [marginAxes marginAxes -2*marginAxes -2*marginAxes] ;
                img0 = zeros(vidRes(2),vidRes(1),3) ;
                imgPreview = image(zeros(vidRes(2),vidRes(1),3)) ;
            % Axes formatting
                axis equal
                axis off
                %zoom on
                axPreview.YDir = 'reverse' ;
                axPreview.XTick = [] ;
                axPreview.YTick = [] ;
                axPreview.XLim = [0 vidRes(1)] ;
                axPreview.YLim = [0 vidRes(2)] ;
        % Add a zoom button
            zoomBtn = uicontrol('style','togglebutton',...
                                'string','zoom',...
                                'units','normalized',...
                                'position',[.95 .955 .035 .025],...
                                'foregroundcolor','r',...
                                'fontweight','bold') ;
            zoomBtn.Callback = @(src,evt)zoomCallBack(src) ;
        % Add a DERIVATIVE btn
            PREVIEW.derivBtn = uicontrol('style','popupmenu',...
                                'units','normalized',...
                                'position',[.935 .919 .05 .025],...
                                'value',1) ;
            PREVIEW.derivBtn.String = {'flat','gradient','laplacian','noise level'} ;
        % Setup imRectangle for ROI
            rectROI = imrect(axPreview,INFOS.Custom.ROI.Value) ; %[0 0 vidRes+1]
            set(rectROI,'Tag','rectROI') ;
            fcn = makeConstrainToRectFcn('imrect',axPreview.XLim, axPreview.YLim) ;
            setPositionConstraintFcn(rectROI,fcn) ;
            addNewPositionCallback(rectROI,@(pos)updateSettings(rectROI,pos)) ;
            rectROI.Deletable = false ;
        % Setup the histogram
            axHisto = axes() ;
                axHisto.Position = [marginMenu marginMenu menuSize/(1+menuSize)-2*marginMenu histoHeight] ;
                axHisto.XTick = [] ;
                axHisto.YTick = [] ;
                box on ;
        % Create the menu
            createMenu() ;
        % Init UI Values ;
            initUIValues() ;
        % Launch the preview
            PREVIEW.LastFrame = img0 ;
            PREVIEW.AxHistObj = axHisto ;
            PREVIEW.ImgObj = preview(cam.VidObj,imgPreview) ;
            setappdata(PREVIEW.ImgObj,'UpdatePreviewWindowFcn',@(obj,event,hImage)updatePreview(obj,event,hImage)) ;
    end

% ZOOM ON PREVIEW /!\ rect ROI is not settable if zoom is active
    function zoomCallBack(src) 
        if src.Value
            zoom(fig,'on') ;
        else
            zoom(fig,'off') ;
        end
    end


% CLOSE THE FIGURE
    function closeFigure()
        stoppreview(cam.VidObj) ;
        fig.CloseRequestFcn = @(src,evt)closereq ;
        close(fig) ;
        cam.CurrentState = 'connected' ;
    end


        
end