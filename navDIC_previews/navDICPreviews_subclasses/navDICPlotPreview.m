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
                    plotMachin = 'force_strain3D' ;
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
                            prev.YDataSources{1} = ['meanNoNaN((hd.Seeds(end).Displacements(:,','2',',:)),1)'] ;
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
                            prev.XDataSources{1} = ['meanNoNaN((hd.Seeds(1).Strains(:,2,:)),1)'] ;
                            prev.YDataSources{1} = ['meanNoNaN((hd.Seeds(1).Strains(:,1,:)),1)'] ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Poisson') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                        case 'FORCE_STRAIN'
                            prev.XDataSources{1} = ['meanNoNaN((hd.Seeds(end).Strains(:,',num2str(1),',:)),1)'] ;
                            prev.YDataSources{1} = ['hd.InputData-hd.InputData(1)'] ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Poisson') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                        case 'FORCE_STRAIN3D'
                            nb = length(hd.InputData-hd.InputData(1)) ;
                            for i = 1:length(hd.Seeds(end))
                                prev.XDataSources{1+i-1} = ['reshape(hd.Seeds(end).Strains(1,1,:,1),[',num2str(nb),' 1])'] ;
                                prev.YDataSources{1} = 'hd.InputData-hd.InputData(1)' ;
                                prev.XDataSources{2} = ['reshape(hd.Seeds(end).Strains(1,1,:,2),[',num2str(nb),' 1])'] ;
                                prev.YDataSources{2} = 'hd.InputData-hd.InputData(1)' ;
                                prev.lines(1) = plot(NaN,NaN,'tag','Cam 1') ;
                                prev.lines(2) = plot(NaN,NaN,'tag','Cam 2') ;
                                prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                                prev.timeMarkers(1) = plot(NaN,NaN) ;
                                prev.timeMarkers(2) = plot(NaN,NaN) ;
                            end
                    end
                    title(regexprep(plotMachin,{'_'},{' '}))
                    prev = updatePreview(prev,hd) ;
                    if strcmpi(plotMachin,'FORCE_STRAIN3D')
                        set(prev.lines(1),'linestyle','-.','linewidth',1,'marker','.','markersize',20,'color','b')
                        set(prev.lines(2),'linestyle','-.','linewidth',1,'marker','.','markersize',20,'color','r')
                        set(prev.timeMarkers,'linestyle','none','marker','o','markersize',12,'linewidth',2)
                    else
                        set(prev.lines,'linestyle','-.','linewidth',1,'marker','.','markersize',20,'color','b')
                        set(prev.timeMarkers,'linestyle','none','marker','o','markersize',12,'linewidth',2)
                    end
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
                                if ~isempty(xdata) && ~isempty(ydata)
                                    prev.lines(l).XData = xdata' ;
                                    prev.lines(l).YData = ydata' ;
                                    prev.timeMarkers(l).XData = xdata(hd.CurrentFrame) ;
                                    prev.timeMarkers(l).YData = ydata(hd.CurrentFrame) ;
                                end
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