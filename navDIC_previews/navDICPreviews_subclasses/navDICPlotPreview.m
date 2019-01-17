classdef navDICPlotPreview < navDICPreview
    
    properties
        Axes = [] ;
        lines = gobjects(0) ;
        timeMarkers = gobjects(0) ;
        XDataSources = {} ;
        YDataSources = {} ;
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
                % TEMPORARY CODE =======================
                    % Force vs Time
                    plotMachin = 'poisson' 
                    timeString = 'sum(bsxfun(@times,bsxfun(@minus,hd.TimeLine,hd.TimeLine(1,:)),[0 0 0 3600 60 1]),2)' ;
                    switch upper(plotMachin)
                        case 'FORCE_TIME'
                            prev.XDataSources{1} = timeString ;
                            prev.YDataSources{1} = 'hd.InputData-hd.InputData(1)' ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Force/Time') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                        case 'DISP_TIME'
                            prev.XDataSources{1} = timeString ;
                            prev.YDataSources{1} = 'meanNoNaN((hd.Seeds(1).Displacements(:,2,:)),1)' ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Displacement/Time') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                        case 'STRAIN_TIME'
                            for e = 1:3
                                prev.XDataSources{e} = timeString ;
                                prev.YDataSources{e} = ['meanNoNaN((hd.Seeds(end).Strains(:,',num2str(e),',:)),1)'] ;
                                prev.lines(e) = plot(NaN,NaN,'tag','Srain/Time') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                                prev.timeMarkers(e) = plot(NaN,NaN) ;
                            end
                        case 'POSITION_TIME'
                            prev.XDataSources{1} = timeString ;
                            prev.YDataSources{1} = ['meanNoNaN((hd.Seeds(end).Displacements(:,','1',',:)),1)'] ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Position/Time') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                        case 'VELOCITY_TIME'
                            prev.XDataSources{1} = timeString ;
                            prev.YDataSources{1} = ['cat(3,0,diff(meanNoNaN((hd.Seeds(end).Displacements(:,','1',',:)),1)))'] ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Position/Time') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                        case 'POISSON'
                            prev.XDataSources{1} = ['meanNoNaN((hd.Seeds(end).Strains(:,2,:)),1)'] ;
                            prev.YDataSources{1} = ['meanNoNaN((hd.Seeds(end).Strains(:,1,:)),1)'] ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Poisson') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                        case 'FORCE_STRAIN'
                            prev.XDataSources{1} = ['meanNoNaN((hd.Seeds(1).Strains(:,2,:)),1)'] ;
                            prev.YDataSources{1} = ['hd.InputData-hd.InputData(1)'] ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Poisson') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                    end
                    title(regexprep(plotMachin,{'_'},{' '}))
                    prev = updatePreview(prev,hd) ;
                    set(prev.lines,'linestyle','-.','linewidth',1,'marker','.','markersize',20)
                    set(prev.timeMarkers,'linestyle','none','marker','o','markersize',12,'linewidth',4)
                %=======================================
            end
            
        % UPDATE
            function prev = updatePreview(prev,hd)
                % Superclass updating function
                    prev = updatePreview@navDICPreview(prev,hd) ;
                    if ~prev.isValid ; return ; end
                    % TEMPORARY CODE =======================
                        for l = 1:length(prev.lines)
                            xdata = [] ;
                            ydata = [] ;
                            try
                                xdata = eval(prev.XDataSources{l}) ;
                                ydata = eval(prev.YDataSources{l}) ;
                            end
                            if ~isempty(xdata) && ~isempty(ydata)
                                prev.lines(l).XData = xdata(1:length(ydata)) ;
                                prev.lines(l).YData = ydata ;
                                prev.timeMarkers(l).XData = xdata(hd.CurrentFrame) ;
                                prev.timeMarkers(l).YData = ydata(hd.CurrentFrame) ;
                            end
                        end
                    %=======================================
            end
        
        % DESTRUCTOR
            function closePreview(obj)
            end
            
    end
    
            
    % OTHER FUNCTIONS
    
        methods(Static)                
        end
    
end