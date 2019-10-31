function H = drawingTool(varargin)
% varargin is used in cases where no base rectangle is needed.
% varargin parameters
%    'background', image
%    'baseframe', boolean (draw a baseframe)
%    'allowed shapes', cell array of strings
%           - rect+
%           - rect-
%           - elli+
%           - elli-
%           - poly+
%           - poly-
%           - point
%           - line
%           - polyline
%           - spline
%    'drawROI', boolean
%    'updateCallback', function @(handles)fun(handles)
%    'uiContextMenu', {parameter,{allowed values},default value}


    %global handles (for debugging) 
    
    % INPUT INIT
        parseInputs(varargin{:}) ;
        H ;
        
    % FIGURE INITIALIZATION  
        initFig() ;
        updateDrawing(0) ;
        
    % WAIT FOR THE USER TO CLOSE THE FIGURE
        while isvalid(H.Figure)
            drawnow ;
            %pause(0.04) ;
        end

    % RETURN HANDLES
        
% ===================================================================================================================    
% CHILD FUNCTIONS
% ===================================================================================================================


% ADD A SHAPE
    function newShape(shape,msg,varargin)
        numShape = length(H.Shapes)+1 ;
        if regexp(msg,'impolyline','once') ; shape = 'impolyline' ; end
        %msg = [msg , ' [',shape,']'] ;
        % Shape-depending properties
            switch shape
                case 'imrect'
                    SHAPE = imrect(H.Axes,varargin{:}) ; 
                case 'imellipse'
                    SHAPE = imellipse(H.Axes,varargin{:}) ;
                case 'impoly'
                    SHAPE = impoly(H.Axes,varargin{:}) ;
                case 'impoint'
                    SHAPE = impoint(H.Axes,varargin{:}) ;
                case 'imline'
                    SHAPE = imline(H.Axes,varargin{:}) ;
                case 'impolyline'
                    SHAPE = impoly(H.Axes,varargin{:},'Closed',false) ;
                    shape = 'impoly' ;
                    msg = 'impolyline' ;
                otherwise
                    return ;
            end
            if isempty(SHAPE) ; return ; end
            H.Shapes{numShape} = SHAPE ;
            fcn = makeConstrainToRectFcn(shape,H.Axes.XLim, H.Axes.YLim) ;
        % Type-depending properties
            switch msg(1)
                case '+'
                    setColor(H.Shapes{numShape},'b') ;
                case '-'
                    setColor(H.Shapes{numShape},'r') ;
            end
        % Common properties 
            set(H.Shapes{numShape},'DisplayName',msg) ;
            setPositionConstraintFcn(H.Shapes{numShape},fcn) ;
            addNewPositionCallback(H.Shapes{numShape},@(pos)updateDrawing(numShape)) ;
        % /!\ this listener has to be cleared before figure closing
            H.Listeners{numShape} = addlistener(H.Shapes{numShape},'ObjectBeingDestroyed',@(src,evt)updateDrawing(0)) ;
        % UIContextMenu
            children = get(H.Shapes{numShape},'children') ;
            contextmenu = get(children(1),'uicontextmenu') ;
            if isempty(contextmenu) ; contextmenu = get(H.Shapes{numShape},'uicontextmenu') ; end 
            menu = uimenu(contextmenu,'separator','on','label','No Custom Data','enable','off') ;
            % If a context menu has been asked for
                if ~isempty(H.uiContextMenu)
                    menu.Label = 'Custom Data' ;
                    menu.Enable = 'on' ;
                    for m = 1:size(H.uiContextMenu,1)
                        submenu = uimenu(menu) ;
                        submenu.Label = H.uiContextMenu{m,1} ;
                        submenu.Callback = @(src,evt)uiContextCallback(src,evt,numShape,'menu') ;
                        for sm = 1:length(H.uiContextMenu{m,2})
                            subsubmenu = uimenu(submenu) ;
                            subsubmenu.Label = H.uiContextMenu{m,2}{sm} ;
                            subsubmenu.Callback = @(src,evt)uiContextCallback(src,evt,numShape,'submenu') ;
                            if strcmp(strtrim(H.uiContextMenu{m,2}{sm}),strtrim(H.uiContextMenu{m,3}))
                                subsubmenu.Checked = 'on' ;
                                if H.initFinished
                                    H.uiContextMenuData(numShape).(submenu.Label) = strtrim(subsubmenu.Label) ;
                                end
                            end
                        end
                    end
                end
        % Update
            if H.initFinished
                updateDrawing(0) ;
            end
        % Depending on the shape, loop
            if isempty(varargin)
                switch shape
                    case 'impoint'
                        newShape(shape,msg) ;
                    otherwise
                end
            end
    end



% UPDATE THE DRAWING
    function updateDrawing(numShape)
        % Save the shape that called the function
            if numShape>0
                H.numShape = numShape ;
            end
        % Update ROI if needed
            if H.drawROI
                computeROI() ;
                H.SrfROI.AlphaData = (1-H.ROI)*0.1 ;
            end
        % Backup Geometries
            H.Geometries = [] ;
            for s = 1:length(H.Shapes)
                if ~isvalid(H.Shapes{s}) 
                    H.Geometries(end+1).isValid = false ;
                    continue ; 
                end
                H.Geometries(end+1).Class = class(H.Shapes{s}) ;
                H.Geometries(end).Bool = get(H.Shapes{s},'DisplayName') ;
                H.Geometries(end).Position = getPosition(H.Shapes{s}) ;
                H.Geometries(end).isValid = true ;
            end
        % Callback
            if strcmp(H.btnCallback.State,'on') 
                executeCallbackFcn() ;
            end
    end


% COMPUTE THE ROI
    function computeROI()
        % Basically computes boolean operations
            ROI = logical(H.BackGround(:,:,1).*0) ;
            shapes = H.Shapes ;
            for s = 1:length(shapes)
                try
                    msg = get(shapes{s},'DisplayName') ;
                    switch msg
                        case '+'
                            ROI = ROI | logical(createMask(shapes{s},H.SrfROI)) ;
                        case '-'
                            ROI = ROI & logical((1-createMask(shapes{s},H.SrfROI))) ;
                        otherwise
                    end
                end
            end
            H.ROI = double(ROI) ;
    end

% EXECUTE THE USER-SPECIFIED CALLBACK FUNCTION
    function executeCallbackFcn() 
        H = H.updateCallback(H) ;
    end

% CALLABACK FOR OPTIONAL UIMENUS
    function uiContextCallback(src,evt,numShape,menuType)
            switch menuType
                case 'menu'
                    % Update checks
                        child = src.Children ;
                        Value = H.uiContextMenuData(numShape).(src.Label) ;
                        isChecked = ismember(strtrim({child.Label}),Value) ;
                        set(child(~isChecked),'Checked','off') ;
                        set(child(isChecked),'Checked','on') ;
                case 'submenu'
                    % Get infos
                        Value = strtrim(src.Label) ;
                        Label = strtrim(src.Parent.Label) ;
                    % Set Handles Data
                        H.uiContextMenuData(numShape).(Label) = Value ;
                    % Change the checked status
                        child = src.Parent.Children ;
                        set(child,'Checked','off') ;
                        set(src,'checked','on') ;
                    % Execute the callback function
                        updateDrawing(numShape) ;
            end
    end


% EXECUTE IF CLICK ON FIGURE
    function clickOnFigure(src,evt)
        return ;
        imgPt = H.Axes.CurrentPoint(1,1:2) ;
        switch H.Figure.SelectionType
            case 'normal'
                disp('left click')
                return ;
            case 'alt'
                disp('right click')
            case 'open'
                disp('double click')
            otherwise
                disp('not known button')
        end
    end

         
% INITIALIZATION OF THE MAIN FIGURE
    function initFig()
        
        marg = 3 ; % Pixels
        refImg = H.BackGround ;
        figRelSize = .55 ;

        % Figure
            H.Figure = figure('tag','drawingToolFig'...
                                ,'Name','drawingTool'...
                                ,'NumberTitle','off') ;
            H.Figure.CloseRequestFcn = @(src,evt)closeFigure() ;
            % Fix the aspect ratio to image
                maxSize = get(0,'defaultfigureposition') ; maxSize = maxSize(3:4)*figRelSize  ;
                Ratio = max((size(refImg')+2*marg)./(maxSize-[0 55]));
                H.Figure.Position([4,3]) = round(size(refImg)/Ratio)+2*marg ;
            % Center the figure on screen
                H.Figure.Units = 'normalized' ;
                H.Figure.Position(1:2) = [.5-H.Figure.Position(3:4)/2] ;

        % Image    
            H.Axes = axes('position',[0 0 1 1]) ;
                H.Axes.Units = 'pixels' ;
                H.Axes.Position = H.Axes.Position + [marg marg -2*marg -2*marg] ;
                H.Axes.Units = 'normalized' ;
                H.Axes.YDir = 'reverse' ;
                H.Img = imagesc(refImg) ;
                H.Img.CData = repmat(refImg,[1 1 3]) ;
                axis equal
                axis tight
                axis off

        % ToolBar    
            H.ToolBar = uitoolbar(H.Figure) ;
            addButtons ;  
            
        % Click callback
            addlistener(H.Figure,'WindowMousePress',@(src,evt)clickOnFigure(src,evt)) ;
                
        % ROI
            if H.drawROI
                H.ROI = H.BackGround(:,:,1).*0+1 ;
                H.SrfROI = imagesc(refImg) ;
                H.SrfROI.CData = repmat(H.ROI,[1 1 3]).*0+255 ;
            end

        % Basis Rectangle if needed
            if H.BaseFrame
                newShape('imrect','+',[0 0 size(H.Img.CData,2) size(H.Img.CData,1)]) ;
                H.Shapes{end}.Deletable = 0 ;
            end
            
        % Create geometries that are existing (handle struct was given as  input)
            for g = 1:length(H.Geometries)
                if ~H.Geometries(g).isValid ; continue ; end
                newShape(H.Geometries(g).Class,...
                            H.Geometries(g).Bool,...
                            H.Geometries(g).Position) ;
            end
            
        % Go !
            H.initFinished = true ;
        
    end



% ADD THE BUTTONS
    function addButtons
        
        % Buttons as cell array
            buttons = H.AllowedShapes ;
            if ~isempty(H.updateCallback) ; buttons{end+1} = 'update' ; end
            buttons{end+1} = 'OK' ;
            
        % Button declaration    
            for btn = 1:length(buttons) 
                % Initialize the Button
                    theBtn = uipushtool(H.ToolBar) ;
                    icon = [] ;
                    switch buttons{btn}
                        case 'rect+'
                            theBtn.TooltipString = 'Add a Rectangle' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('imrect','+',[]) ;
                            icon = 'iconAddRect.png' ;
                        case 'rect-'
                            theBtn.TooltipString = 'Remove a Rectangle' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('imrect','-',[]) ;
                            icon = 'iconRemRect.png' ;
                        case 'elli+'
                            theBtn.TooltipString = 'Add an Ellipse' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('imellipse','+',[]) ;
                            icon = 'iconAddElli.png' ;
                        case 'elli-'
                            theBtn.TooltipString = 'Remove an Ellipse' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('imellipse','-',[]) ;
                            icon = 'iconRemElli.png' ;
                        case 'poly+'
                            theBtn.TooltipString = 'Add a Polygon' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('impoly','+',[]) ;
                            icon = 'iconAddPoly.png' ;
                        case 'poly-'
                            theBtn.TooltipString = 'Remove a Polygon' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('impoly','-',[]) ;
                            icon = 'iconRemPoly.png' ;
                        case 'point'
                            theBtn.TooltipString = 'Draw a Point' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('impoint','+',[]) ;
                        case 'line'
                            theBtn.TooltipString = 'Draw a Line' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('imline','+',[]) ;
                        case 'polyline'
                            theBtn.TooltipString = 'Draw a Polyline' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('impolyline','+',[]) ;
                        case 'spline'
                            theBtn.TooltipString = 'Draw a Spline' ;
                            theBtn.ClickedCallback = @(src,evt)newShape('Spline','+',[]) ;
                        case 'update'
                            % ToggleTool
                                delete(theBtn) ;
                                theBtn = uitoggletool(H.ToolBar) ;
                                H.btnCallback = theBtn ;
                                theBtn.State = 'off' ;
                            theBtn.TooltipString = 'Auto Update' ;
                            theBtn.ClickedCallback = @(src,evt)updateDrawing(0) ;
                            icon = 'help_fx.png' ;
                        case 'OK'
                            theBtn.TooltipString = 'Validate the Draw and Quit' ;
                            theBtn.ClickedCallback = @(src,evt)closeFigure() ;
                            icon = 'iconRoiOK.png' ;
                    end
                % Draw The Icon
                    if isempty(icon) ; icon = 'tool_shape_stroke.png' ; end
                    [icon,~,alpha] = imread(icon) ; 
                    icon = im2single(icon) ;
                    icon(repmat(alpha,[1 1 3])==0) = NaN ;
                    theBtn.CData = icon ;
            end

    end

% CLOSE THE FIGURE
    function closeFigure()
            H.Figure.CloseRequestFcn = @(src,evt)closereq ;
        % Delete Listeners
            try
                for l = 1:length(H.Listeners)
                    H.Listeners{l}.Enabled = false ;
                end
            end
        % Close figure
            close(H.Figure) ;
    end

% INIT HANDLES WITH INPUTS
    function parseInputs(varargin)
        % Was a handle structure given ?
            hasToBeInit = true ;
            if ~isempty(varargin) && strcmp(class(varargin{1}),'struct')
                H = varargin{1} ;
                hasToBeInit = false ;
            end
        % Init if needed
            if hasToBeInit 
                % Initialize Handles
                    H = {} ;
                    H.Type = 'drawingToolHandles' ;
                    H.BackGround = ones(480,640)*0.5 ;
                    H.BaseFrame = false ;
                    H.AllowedShapes = {'rect+','rect-','elli+','elli-','poly+','poly-','point','line','polyline','spline'} ;
                    H.drawROI = false ;
                    H.updateCallback = @(H)H ;
                    H.uiContextMenu = [] ;
                    H.uiContextMenuData = [] ;
                % VARARGIN PROCESS
                    for arg = 1:2:length(varargin)
                        switch upper(varargin{arg})
                            case 'BACKGROUND'
                                H.BackGround = varargin{arg+1} ;
                            case 'BASEFRAME'
                                H.BaseFrame = varargin{arg+1} ;
                            case 'ALLOWEDSHAPES'
                                H.AllowedShapes = varargin{arg+1} ;
                            case 'DRAWROI'
                                H.drawROI = varargin{arg+1} ;
                            case 'UPDATECALLBACK'
                                H.updateCallback = varargin{arg+1} ;
                                H.CallbackResult = [] ;
                            case 'UICONTEXTMENU'
                                H.uiContextMenu = varargin{arg+1} ;
                        end
                    end
                H.Geometries = [] ;
                H.CallbackResult = [] ;
            end
        % HANDLES BASE SETTING
            H.numShape = 0 ;
            H.Shapes = {} ;
            H.Listeners = {} ;
            H.initFinished = false ;
    end


end
