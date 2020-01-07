classdef navDICPlotPreview < navDICPreview
    
    properties
        Axes = gobjects(0)
        Curves = gobjects(0)
        TimeMarkers = gobjects(0)
        XDataSources = {}
        YDataSources = {}
    end
    
    properties (Hidden) % UI controls
        CurvePanel = gobjects(0)
        AddCurveBtn = gobjects(0)
        RemoveCurveBtn = gobjects(0)
        DupplicateCurveBtn = gobjects(0)
        CurveList = gobjects(0)
        CurveXData = gobjects(0)
        CurveYData = gobjects(0)
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
                    else
                        prev.CurveXData.String = prev.XDataSources{crv} ;
                        prev.CurveYData.String = prev.YDataSources{crv} ;
                    end
                % Update curves
                    for l = 1:length(prev.Curves)
                        [xdata,ydata] = str2data(prev,hd,prev.XDataSources{l},prev.YDataSources{l}) ;
                        if ~isempty(xdata) && ~isempty(ydata)
                            % Apply to curves
                                prev.Curves(l).XData = xdata(1:length(ydata)) ;
                                prev.Curves(l).YData = ydata ;
                                currentFrame = min(numel(ydata),hd.CurrentFrame) ;
                                prev.TimeMarkers(l).XData = xdata(currentFrame) ;
                                prev.TimeMarkers(l).YData = ydata(currentFrame) ;
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
                        seedStr = ['hd.Seeds(',num2str(seedID),')'] ;
                    end
                % Prepare the data
                    timeString = 'sum(bsxfun(@times,bsxfun(@minus,hd.TimeLine,hd.TimeLine(1,:)),[0 0 0 3600 60 1]),2)' ;
                    switch upper(preset)
                        case 'INPUT VS. TIME'
                            prev.addCurve('Input/Time' ... % Name
                                            ,'hd.InputData-hd.InputData(1)' ... % YData
                                            , timeString ... % XData
                                            ) ;
                        case 'DISPLACEMENT VS. TIME'
                            for fi = {'u1','u2'}
                                prev.addCurve([fi{:}],[seedStr '.DataFields.' [fi{:}]],timeString) ;
                            end
                        case 'STRAIN VS. TIME'
                            for fi = {'L11','L22','L12'}
                                prev.addCurve([fi{:}],[seedStr '.DataFields.' [fi{:}]],timeString) ;
                            end
                        case 'POSITION VS. TIME'
                            for fi = {'x1','x2'}
                                prev.addCurve([fi{:}],[seedStr '.DataFields.' [fi{:}]],timeString) ;
                            end
                        case 'VELOCITY VS. TIME'
                            for fi = {'v1','v2'}
                                prev.addCurve([fi{:}],[seedStr '.DataFields.' [fi{:}]],timeString) ;
                            end
                        case 'POISSON'
                            prev.addCurve('Poisson' ...
                                            ,['-' seedStr '.DataFields.TSe2'] ...
                                            ,[seedStr '.DataFields.TSe1']) ;
                        case 'INPUT VS. STRAIN'
                            prev.addCurve('Input/Strain' ...
                                            ,'hd.InputData' ...
                                            ,[seedStr '.DataFields.Le1']) ;
                        case 'TRUE STRESS VS. TRUE STRAIN'
                            prev.addCurve('TrueStress/TrueStrain' ...
                                            ,['(1/(1*1)).*reshape(hd.InputData,1,1,[]).*' seedStr '.DataFields.lambda1' ] ...
                                            ,[seedStr '.DataFields.TSe1']) ;
                        case 'CUSTOM'
                            prev.addCurve(['curve.' num2str(numel(prev.CurveList.String)+1)]) ; 
                    end
                    %title(regexprep(preset,{'_'},{' '}))
                    prev = updatePreview(prev,hd) ;
            end
        
        % TRANSLATE STRING TO DATA
            function [X,Y] = str2data(obj,hd,Xstr,Ystr)
                X = [] ; Y = [] ;
                % Process the strings with keywords and replacement
                    keywords = {} ; replace = {} ;
                    % Time line
                        keywords{end+1} = '@time' ; replace{end+1} = 'sum(bsxfun(@times,bsxfun(@minus,hd.TimeLine,hd.TimeLine(1,:)),[0 0 0 3600 60 1]),2)' ; 
                        for s = 1:numel(hd.Seeds)
                            keywords{end+1} = ['@' hd.Seeds(s).Name '.'] ; 
                            replace{end+1} = ['hd.Seeds(' num2str(s) ').DataFields.'] ;
                        end
                    % Process
                        Xstr = regexprep(Xstr,keywords,replace) ;
                        Ystr = regexprep(Ystr,keywords,replace) ;
                % Try retrieving the data
                    try
                        X = eval(Xstr) ;
                        Y = eval(Ystr) ;
                    catch error
                        warning(error.message) ;
                    end
                % Format the data
                    if ~isempty(X) && ~isempty(Y)
                        % To 2D
                            X = X(:,:) ;
                            Y = Y(:,:) ;
                        % To 1D ;
                            % Column vector to row-vector
                                if size(X,2)==1 ; X = X(:)' ; end
                                if size(Y,2)==1 ; Y = Y(:)' ; end
                            % Mean over the 1st dim.
                                X = mean(X,1,'omitnan') ;
                                Y = mean(Y,1,'omitnan') ;
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
                prev.CurveXData = uicontrol(prev.CurvePanel ...
                                            ,'style','edit' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[btnWidth+listWidth .5 1-btnWidth-listWidth .5]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)changeDataSource(prev) ...
                                        ) ; 
                prev.CurveYData = uicontrol(prev.CurvePanel ...
                                            ,'style','edit' ...
                                            ,'units','normalized' ...
                                            ,'OuterPosition',[btnWidth+listWidth 0 1-btnWidth-listWidth .5]  + [1 1 -2 -2]*margin ...
                                            ,'Callback',@(src,evt)changeDataSource(prev) ...
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
            function addCurve(prev,name,ydata,xdata,update)
                if nargin<2 ; prev.addPreset ; return ; end
                if nargin<3 ; ydata = 'NaN(hd.nFrames,1)' ; end
                if nargin<4 ; xdata = '(1:hd.nFrames)''' ; end
                if nargin<5 ; update = true ; end
                % Add the curve data
                    prev.CurveList.String{end+1} = name ;
                    prev.XDataSources{end+1} = xdata ;
                    prev.YDataSources{end+1} = ydata ;
                    prev.CurveList.Value = numel(prev.CurveList.String) ;
                % Init the curve plots
                    clrOrder = prev.Axes.ColorOrderIndex ;
                    prev.Curves(end+1) = plot(prev.Axes ...
                                                ,NaN,NaN ...
                                                ,'linestyle','-' ...
                                                ,'linewidth',1 ...
                                                ,'marker','.' ...
                                                ,'markersize',7 ...
                                                ,'tag',name ...
                                                ,'displayname',name ...
                                                ) ;
                    prev.Axes.ColorOrderIndex = clrOrder ;
                    prev.TimeMarkers(end+1) = plot(prev.Axes ...
                                                    ,NaN,NaN ...
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
                % Add the new curve
                    addCurve(prev,name,ydata,xdata,update)
            end
            function removeCurve(prev,crv,update)
                if nargin<2 ; crv = prev.CurveList.Value ; end
                if nargin<3 ; update = true ; end
                % Remove data sources
                    prev.XDataSources(crv) = [] ;
                    prev.YDataSources(crv) = [] ;
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
                prev.updatePreview ;
            end
    end
            
    % OTHER FUNCTIONS
    
        methods(Static)                
        end
    
end