classdef navDICPlotPreview < navDICPreview
    
    properties
        Axes = gobjects(0)
        Curves = gobjects(0)
        TimeMarkers = gobjects(0)
        XDataSources = {}
        YDataSources = {}
        ZDataSources = {}
    end
    
    properties (Hidden) % UI controls
        CurvePanel = gobjects(0)
        AddCurveBtn = gobjects(0)
        RemoveCurveBtn = gobjects(0)
        DupplicateCurveBtn = gobjects(0)
        CurveList = gobjects(0)
        CurveType = gobjects(0)
        CurveXData = gobjects(0)
        CurveYData = gobjects(0)
        CurveZData = gobjects(0)
        Legend = gobjects(0)
    end
    
    methods
        
        % CONSTRUCTOR        
            function prev = navDICPlotPreview(hd,varargin)
                % Superclass constructor call
                    prev = prev@navDICPreview(hd,varargin{:}) ;
                    prev.fig.Name = 'navDIC Plot Preview' ;
                % Reset the toolbar
                    prev.fig.ToolBar = 'figure' ;
                    prev.fig.MenuBar = 'figure' ;
                % Set the axes
                    prev.Axes = axes('outerposition',[0 0 1 1]) ;
                        box(prev.Axes,'on')
                        grid(prev.Axes,'on')
                        xlabel(prev.Axes,'XData') ;
                        ylabel(prev.Axes,'YData') ;
                        zlabel(prev.Axes,'ZData') ;
                    prev.Legend = legend(prev.Axes) ;
                        prev.Legend.EdgeColor = 'k' ;
                        prev.Legend.Color = 'none' ;
                        prev.Legend.Location = 'best' ;
                % Init the curve panel
                    prev.initCurvePanel ;
                    clearPreview(prev) ;
                % Add a curve
                    prev.addPreset ;
            end
            
        % UPDATE
            function prev = updatePreview(prev,hd)
                if nargin<2 ; global hd ; end
                % Superclass updating function
                    prev = updatePreview@navDICPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
                % Update the panel
                    crv = prev.CurveList.Value ;
                    if isempty(prev.XDataSources)
                        prev.CurveXData.String = '' ;
                        prev.CurveYData.String = '' ;
                        prev.CurveZData.String = '' ;
                    else
                        prev.CurveXData.String = prev.XDataSources{crv} ;
                        prev.CurveYData.String = prev.YDataSources{crv} ;
                        prev.CurveZData.String = prev.ZDataSources{crv} ;
                    end
                % Update curves
                    for l = 1:length(prev.Curves)
                        [xdata,ydata,zdata] = str2data(prev,hd,prev.XDataSources{l},prev.YDataSources{l},prev.ZDataSources{l}) ;
                        if ~isempty(xdata) && ~isempty(ydata) && ~isempty(zdata)
                            % Format the data
                                currentFrame = min(size(ydata,1),hd.CurrentFrame) ;
                                xdata = xdata + 0*ydata + 0*zdata ;
                                ydata = ydata + 0*xdata + 0*zdata ;
                                zdata = zdata + 0*xdata + 0*ydata ;
                                % Process
                                    if size(xdata,2)>1 
                                            switch prev.CurveType.String{prev.CurveType.Value}
                                                case 'all' % Do nothing
                                                case 'mean'
                                                    xdata = mean(xdata,2,'omitnan') ;
                                                    ydata = mean(ydata,2,'omitnan') ;
                                                    zdata = mean(zdata,2,'omitnan') ;
                                                case 'median'
                                                    xdata = median(xdata,2,'omitnan') ;
                                                    ydata = median(ydata,2,'omitnan') ;
                                                    zdata = median(zdata,2,'omitnan') ;
                                                case 'confidence'
                                                    xdata = mean(xdata,2,'omitnan')+[-1 1].*std(xdata,1,2,'omitnan') ;
                                                    ydata = mean(ydata,2,'omitnan')+[-1 1].*std(ydata,1,2,'omitnan') ;
                                                    zdata = mean(zdata,2,'omitnan')+[-1 1].*std(zdata,1,2,'omitnan') ;
                                                case 'minmax'
                                                    xdata = [min(xdata,[],2) max(xdata,[],2)] ;
                                                    ydata = [min(ydata,[],2) max(ydata,[],2)] ;
                                                    zdata = [min(zdata,[],2) max(zdata,[],2)] ;
                                                otherwise
                                            end
                                    end
                                % Concatenate
                                    nans = NaN(size(ydata(1,:))) ;
                                    xdata = cat(1,xdata,nans) ; ydata = cat(1,ydata,nans) ; zdata = cat(1,zdata,nans) ;
                            % Apply to curves
                                prev.Curves(l).XData = xdata(:) ;
                                prev.Curves(l).YData = ydata(:) ;
                                prev.Curves(l).ZData = zdata(:) ;
                                prev.TimeMarkers(l).XData = xdata(currentFrame,:) ;
                                prev.TimeMarkers(l).YData = ydata(currentFrame,:) ;
                                prev.TimeMarkers(l).ZData = zdata(currentFrame,:) ;
                        end
                    end
            end
        
        % DESTRUCTOR
            function closePreview(prev)
            end
            
        % PRESET DEFINITION
            function addPreset(prev,preset)
                global hd
                % Let the user choose the plot
                    if nargin<2
                        availablePlots = {'Custom' ...
                                            ,'Position vs. Time' ...
                                            ,'Displacement vs. Time' ...
                                            ,'Velocity vs. Time' ...
                                            ,'Strain vs. Time' ...
                                            ,'Strain Tensor' ...
                                            ,'Poisson' ...
                                            ,'Input vs. Time' ...
                                            ,'Input vs. Strain' ...
                                            ,'True Stress vs. True Strain' ...
                                            } ;
                        [IDs,valid] = listdlg('PromptString','Select A Plot Option :',...
                                                    'initialValue',1,...
                                                    'ListString',availablePlots) ;
                        if ~valid ; IDs = 1 ; end
                        preset = availablePlots{IDs} ;
                    end
                % If needed, choose a seed
                    if ~ismember({'INPUT VS. TIME','CUSTOM'},upper(preset))
                        [seedID,valid] = selectSeeds(hd,'single') ;
                        if ~valid ; return ; end
                        %seedStr = ['hd.Seeds(',num2str(seedID),').DataFields.'] ;
                        seedStr = ['@' hd.Seeds(seedID).Name] ;
                    end
                % Prepare the data
                    timeString = '@time' ;
                    switch upper(preset)
                        case 'INPUT VS. TIME'
                            if isempty(hd.InputData)
                                if isempty(hd.DAQInputs) || isempty(hd.DAQInputs.Inputs) % Do nothing
                                else
                                    for in = 1:numel(hd.DAQInputs.Inputs)
                                        prev.addCurve([hd.DAQInputs.Inputs(in).Name ' vs. Time'] ... % Name
                                                        ,['hd.InputData(:,' num2str(in) ')'''] ... % YData
                                                        , timeString ... % XData
                                                        ) ;
                                    end
                                end
                            elseif istable(hd.InputData) % Table Data
                                for in = 1:size(hd.InputData,2)
                                    varName = hd.InputData.Properties.VariableNames{in} ;
                                    varUnit = '' ;
                                    if ~isempty(hd.InputData.Properties.VariableUnits)
                                        varUnit = [' (' hd.InputData.Properties.VariableUnits{in} ')'] ;
                                    end
                                    prev.addCurve([varName varUnit] ... % Name
                                                    ,['@' varName] ... % YData
                                                    , timeString ... % XData
                                                    ) ;
                                end
                            else % Array Data
                                for in = 1:size(hd.InputData,2)
                                    prev.addCurve(['Input' num2str(in) ' vs. Time'] ... % Name
                                                    ,['hd.InputData(:,' num2str(in) ')'''] ... % YData
                                                    , timeString ... % XData
                                                    ) ;
                                end
                            end
                        case 'DISPLACEMENT VS. TIME'
                            for fi = {'u1','u2'}
                                prev.addCurve([fi{:}],[seedStr '.' [fi{:}]],timeString) ;
                            end
                        case 'STRAIN TENSOR'
                            prev.addCurve('Strain Tensor',[seedStr '.TS22'],[seedStr '.TS11'],[seedStr '.TS12']) ;
                        case 'STRAIN VS. TIME'
                            for fi = {'L11','L22','L12'}
                                prev.addCurve([fi{:}],[seedStr '.' [fi{:}]],timeString) ;
                            end
                        case 'POSITION VS. TIME'
                            for fi = {'x1','x2'}
                                prev.addCurve([fi{:}],[seedStr '.' [fi{:}]],timeString) ;
                            end
                        case 'VELOCITY VS. TIME'
                            for fi = {'v1','v2'}
                                prev.addCurve([fi{:}],[seedStr '.' [fi{:}]],timeString) ;
                            end
                        case 'POISSON'
                            prev.addCurve('Poisson' ...
                                            ,['-' seedStr '.TSe2'] ...
                                            ,[seedStr '.TSe1']) ;
                        case 'INPUT VS. STRAIN'
                            prev.addCurve('Input/Strain' ...
                                            ,'hd.InputData(:,1)' ...
                                            ,[seedStr '.Le1']) ;
                        case 'TRUE STRESS VS. TRUE STRAIN'
                            prev.addCurve('TrueStress/TrueStrain' ...
                                            ,['(1/(1*1)).*reshape(hd.InputData,1,1,[]).*' seedStr '.lambda1' ] ...
                                            ,[seedStr '.TSe1']) ;
                        case 'CUSTOM'
                            prev.addCurve(['curve.' num2str(numel(prev.CurveList.String)+1)]) ; 
                    end
                    %title(regexprep(preset,{'_'},{' '}))
                    prev = updatePreview(prev,hd) ;
            end
        
        % TRANSLATE STRING TO DATA
            function [X,Y,Z] = str2data(obj,hd,Xstr,Ystr,Zstr)
                X = [] ; Y = [] ; Z = [] ;
                % Process the strings with keywords and replacement
                    keywords = {} ; replace = {} ;
                    % Frame Rate
                        keywords{end+1} = '@framerate' ;
                        replace{end+1} = '[0;1./diff(@time)]' ;
                    % Frames
                        keywords{end+1} = '@frame' ;
                        replace{end+1} = '1:hd.nFrames' ;
                    % Time line
                        keywords{end+1} = '@time' ;
                        if isempty(hd.TimeLine)
                            replace{end+1} = '1:hd.nFrames' ;
                        else
                            replace{end+1} = '(hd.TimeLine-hd.TimeLine(1,:))*flip(cumprod([1;60;60;24;1;1]))' ;
                        end
                    % Seed names
                        for s = 1:numel(hd.Seeds)
                            keywords{end+1} = ['@' hd.Seeds(s).Name '\.'] ; 
                            replace{end+1} = ['hd\.Seeds(' num2str(s) ')\.DataFields\.'] ;
                        end
                    % Input names
                        if ~isempty(hd.DAQInputs) && ~isempty(hd.DAQInputs.Inputs)
                            for in = 1:size(hd.DAQInputs.Inputs,2)
                                varName = hd.DAQInputs.Inputs(in).DataName ;
                                keywords{end+1} = ['@' varName] ; 
                                replace{end+1} = ['hd\.InputData(:,' num2str(in) ')'] ;
                            end
                        end
                        if istable(hd.InputData)
                            for in = 1:size(hd.InputData,2)
                                varName = hd.InputData.Properties.VariableNames{in} ;
                                keywords{end+1} = ['@' varName] ; 
                                replace{end+1} = ['hd\.InputData\.' varName] ;
                            end
                        end
                    % Process
                        Xstr = regexprep(Xstr,keywords,replace) ;
                        Ystr = regexprep(Ystr,keywords,replace) ;
                        Zstr = regexprep(Zstr,keywords,replace) ;
                % Try retrieving the data
                    try
                        X = eval(Xstr) ;
                        Y = eval(Ystr) ;
                        Z = eval(Zstr) ;
                    catch error
                        warning(error.message) ;
                    end
                % Format the data to a 2D matrix ([nFrames nCurves])
                    % To double
                        if istable(X) ; X = table2array(X) ; end
                        if istable(Y) ; Y = table2array(Y) ; end
                        if istable(Z) ; Z = table2array(Z) ; end
                    if ~isempty(X) && ~isempty(Y) && ~isempty(Z)
                        % To 2D
                            X = X(:,:)' ;
                            Y = Y(:,:)' ;
                            Z = Z(:,:)' ;
                        % Row vector to column-vector (for InputData)
                            if size(X,1)==1 ; X = X(:) ; end
                            if size(Y,1)==1 ; Y = Y(:) ; end
                            if size(Z,1)==1 ; Z = Z(:) ; end
                    end
            end
            
    end
    
% CURVES HANDLING
    methods
        % INITIALIZE SOURCES PANEL
            function initCurvePanel(prev)
                margin = 0.002 ;
                panelHeight = 0.1 ;
                btnWidth = 0.04 ;
                listWidth = 0.2 ;
                lblWidth = 0.04 ;
                popupWidth = 0.08 ;
                curveTypes = {'mean','median','confidence','minmax','all'} ;
                prev.CurvePanel = uipanel(prev.fig ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[0 0 1 panelHeight] + [1 1 -2 -2]*margin ...
                                            ...,'Callback',@(src,evt)updatePreview(prev) ...
                                        ) ;
                prev.AddCurveBtn = uicontrol(prev.CurvePanel ...
                                            ,'style','pushbutton' ...
                                            ,'string','add' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[0 2/3 btnWidth 1/3]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)addCurve(prev) ...
                                        ) ; 
                prev.RemoveCurveBtn = uicontrol(prev.CurvePanel ...
                                            ,'style','pushbutton' ...
                                            ,'string','rem' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[0 1/3 btnWidth 1/3]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)removeCurve(prev) ...
                                        ) ; 
                prev.DupplicateCurveBtn = uicontrol(prev.CurvePanel ...
                                            ,'style','pushbutton' ...
                                            ,'string','dup' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[0 0 btnWidth 1/3]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)dupplicateCurve(prev) ...
                                        ) ; 
                prev.CurveList = uicontrol(prev.CurvePanel ...
                                            ,'style','listbox' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[btnWidth 0 listWidth 1]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)updatePreview(prev) ...
                                        ) ;   
                uicontrol(prev.CurvePanel ...
                                ,'style','pushbutton' ...
                                ,'String','XData' ...
                                ,'enable','off' ...
                                ,'units','normalized' ...
                                ,'OuterPosition',[btnWidth+listWidth 2/3 lblWidth 1/3]  + [1 1 -2 -2]*margin ...
                            ) ; 
                uicontrol(prev.CurvePanel ...
                                ,'style','pushbutton' ...
                                ,'String','YData' ...
                                ,'enable','off' ...
                                ,'units','normalized' ...
                                ,'OuterPosition',[btnWidth+listWidth 1/3 lblWidth 1/3]  + [1 1 -2 -2]*margin ...
                            ) ; 
                uicontrol(prev.CurvePanel ...
                                ,'style','pushbutton' ...
                                ,'String','ZData' ...
                                ,'enable','off' ...
                                ,'units','normalized' ...
                                ,'OuterPosition',[btnWidth+listWidth 0 lblWidth 1/3]  + [1 1 -2 -2]*margin ...
                            ) ; 
                prev.CurveXData = uicontrol(prev.CurvePanel ...
                                            ,'style','edit' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[btnWidth+listWidth+lblWidth 2/3 1-btnWidth-listWidth-lblWidth-popupWidth 1/3]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)changeDataSource(prev) ...
                                        ) ; 
                prev.CurveYData = uicontrol(prev.CurvePanel ...
                                            ,'style','edit' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[btnWidth+listWidth+lblWidth 1/3 1-btnWidth-listWidth-lblWidth-popupWidth 1/3]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)changeDataSource(prev) ...
                                        ) ; 
                prev.CurveZData = uicontrol(prev.CurvePanel ...
                                            ,'style','edit' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[btnWidth+listWidth+lblWidth 0 1-btnWidth-listWidth-lblWidth-popupWidth 1/3]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)changeDataSource(prev) ...
                                        ) ;
                prev.CurveType = uicontrol(prev.CurvePanel ...
                                            ,'style','listbox' ...
                                            ,'String',curveTypes ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[1-popupWidth 0 popupWidth 1]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)updatePreview(prev) ...
                                        ) ;
                prev.showCurvePanel ;
            end
        
        % SHOW/HIDE SOURCES PANEL
            function showCurvePanel(prev)
                height = prev.CurvePanel.OuterPosition(4) ;
                prev.Axes.OuterPosition(2) = height ;
                prev.Axes.OuterPosition(4) = 1-height ;
                prev.CurvePanel.Visible = 'on' ;
            end
            function hideCurvePanel(prev)
                prev.Axes.OuterPosition(2) = 0 ;
                prev.Axes.OuterPosition(4) = 1 ;
                prev.CurvePanel.Visible = 'off' ;
            end
            
        % ADD/REMOVE CURVES
            function addCurve(prev,name,ydata,xdata,zdata,update)
                if nargin<2 ; prev.addPreset ; return ; end
                if nargin<3 ; ydata = 'NaN(hd.nFrames,1)' ; end
                if nargin<4 ; xdata = '(1:hd.nFrames)''' ; end
                if nargin<5 ; zdata = '0' ; end
                if nargin<6 ; update = true ; end
                % Add the curve data
                    prev.CurveList.String{end+1} = name ;
                    prev.XDataSources{end+1} = xdata ;
                    prev.YDataSources{end+1} = ydata ;
                    prev.ZDataSources{end+1} = zdata ;
                    prev.CurveList.Value = numel(prev.CurveList.String) ;
                % Init the curve plots
                    clrOrder = prev.Axes.ColorOrderIndex ;
                    prev.Curves(end+1) = plot3(prev.Axes ...
                                                ,NaN,NaN,NaN ...
                                                ,'linestyle','-' ...
                                                ,'linewidth',1 ...
                                                ,'marker','.' ...
                                                ,'markersize',7 ...
                                                ,'tag',name ...
                                                ,'displayname',name ...
                                                ) ;
                    prev.Axes.ColorOrderIndex = clrOrder ;
                    prev.TimeMarkers(end+1) = plot3(prev.Axes ...
                                                    ,NaN,NaN,NaN ...
                                                    ,'linestyle','none' ...
                                                    ,'marker','.' ...
                                                    ,'markersize',25 ...
                                                    ,'linewidth',4 ...
                                                    ,'handleVisibility','off' ...
                                                    ) ;
                % Update
                    if update ; prev.updatePreview ; end
            end
            function dupplicateCurve(prev,crv,update)
                if nargin<2 ; crv = prev.CurveList.Value ; end
                if nargin<3 ; update = true ; end
                % Remove data sources
                    name = ['copy-' prev.CurveList.String{crv}] ;
                    xdata = prev.XDataSources{crv} ;
                    ydata = prev.YDataSources{crv} ;
                    zdata = prev.ZDataSources{crv} ;
                % Add the new curve
                    addCurve(prev,name,ydata,xdata,zdata,update)
            end
            function removeCurve(prev,crv,update)
                if nargin<2 ; crv = prev.CurveList.Value ; end
                if nargin<3 ; update = true ; end
                % Remove data sources
                    prev.XDataSources(crv) = [] ;
                    prev.YDataSources(crv) = [] ;
                    prev.ZDataSources(crv) = [] ;
                % Delete graphical objects
                    delete(prev.Curves(crv)) ; prev.Curves(crv) = [] ;
                    delete(prev.TimeMarkers(crv)) ; prev.TimeMarkers(crv) = [] ;
                % Remove from the curve list
                    prev.CurveList.String(crv) = [] ;
                    prev.CurveList.Value = min(prev.CurveList.Value,numel(prev.CurveList.String)) ;
                % Update Preview
                    if update ; prev.updatePreview ; end
            end
        
        % EMPTY THE PLOT
            function clearPreview(prev)
                for crv = numel(prev.CurveList.String):-1:1
                    prev.removeCurve(crv,false) ;
                end
            end
            
        % CHANGE THE DATA SOURCE
            function changeDataSource(prev)
                crv = prev.CurveList.Value ;
                prev.XDataSources{crv} = prev.CurveXData.String ;
                prev.YDataSources{crv} = prev.CurveYData.String ;
                prev.ZDataSources{crv} = prev.CurveZData.String ;
                prev.updatePreview ;
            end
    end
            
    % OTHER FUNCTIONS
    
        methods(Static)                
        end
    
end