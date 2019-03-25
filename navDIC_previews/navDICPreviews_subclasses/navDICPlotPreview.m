classdef navDICPlotPreview < navDICPreview
    
    properties
        Axes = [] ;
        lines = gobjects(0) ;
        timeMarkers = gobjects(0) ;
        XDataSources = {} ;
        YDataSources = {} ;
        sources = {} ;
        legend = gobjects(0) ;
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
                    plotMachin = 'force_strain_Comp-fil' ;
                    timeString = 'sum(bsxfun(@times,bsxfun(@minus,hd.TimeLine,hd.TimeLine(1,:)),[0 0 0 3600 60 1]),2)' ;
                    if ~strcmpi(plotMachin,'force_strain_Comp-fil')
                        prev = listeSource(prev,hd) ;
                        lst = listdlg('PromptString','Select the output you want to plot : ',...
                                'SelectionMode','multiple',...
                                'ListString',prev.sources) ;
                    end
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
                            for i = lst
                                prev.XDataSources{i} = ['meanNoNaN(hd.Seeds(', num2str(i),').Strains(1,1,:),1)'] ;
                                prev.YDataSources{i} = ['hd.InputData-hd.InputData(1)'] ;
                                prev.lines(i) = plot(NaN,NaN,'tag','Poisson') ;
                                prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                                prev.timeMarkers(i) = plot(NaN,NaN) ;
                            end
                        case 'FORCE_STRAIN3D'
                            prev.XDataSources{1} = ['meanNoNaN( hd.Seeds(end).Strains(:,1,:,1),1) '] ;
                            prev.YDataSources{1} = 'hd.InputData-hd.InputData(1)' ;
                            prev.XDataSources{2} = ['meanNoNaN( hd.Seeds(end).Strains(:,1,:,2),1) '] ;
                            prev.YDataSources{2} = 'hd.InputData-hd.InputData(1)' ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Cam 1') ;
                            prev.lines(2) = plot(NaN,NaN,'tag','Cam 2') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                            prev.timeMarkers(2) = plot(NaN,NaN) ;
                        case 'FORCE_STRAIN3DVS2D'
                            ncam = input('Entrer le numero de la camera a verifier: ') ;
                            prev.XDataSources{1} = ['meanNoNaN( hd.Seeds(end).Strains(:,1,:,',num2str(ncam),'),1)'] ;
                            prev.YDataSources{1} = 'hd.InputData-hd.InputData(1)' ;
                            prev.XDataSources{2} = ['meanNoNaN( hd.Seeds(end).Strains(:,2,:,',num2str(ncam),'),1)'] ;
                            prev.YDataSources{2} = 'hd.InputData-hd.InputData(1)' ;
                            prev.lines(1) = plot(NaN,NaN,'tag','3D') ;
                            prev.lines(2) = plot(NaN,NaN,'tag','2D') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.legend = legend(prev.Axes,{'3D','2D'}) ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                            prev.timeMarkers(2) = plot(NaN,NaN) ;
                        case 'FORCE_STRAIN_COMP-FIL'
                            prev = listeSource(prev,hd) ;
                            lstComp = listdlg('PromptString','Select the composite jauge : ',...
                                'SelectionMode','single','ListString',prev.sources) ;
                            lstFil = listdlg('PromptString','Select the fil jauge : ',...
                                'SelectionMode','single','ListString',prev.sources) ;
                            ncam = inputdlg('Entrer le numero de la camera a verifier: ') ;
                            LcompMm = inputdlg('Entrer la longueur du composite (mm) : ') ;
                            LcompPix = ['hd.Cameras(1).Properties.pixObjratio * ', num2str(LcompMm{1})] ;
                            EpsF = ['repmat(meanNoNaN( hd.Seeds(',num2str(lstFil),').Strains(:,1,:,',num2str(ncam{1}),...
                                '),1),size(hd.Seeds(',num2str(lstComp),').Strains(:,1,1,',num2str(ncam{1}),...
                                ') ) ) '] ;
                            EpsApp = ['hd.Seeds(',num2str(lstComp),').Strains(:,1,:,',ncam{1},')'] ;
                            L0 = ['repmat(hd.Seeds(',num2str(lstComp),').L0(:,',ncam{1},...
                                '),[1 1 size(hd.Seeds(',num2str(lstComp),').Strains,3) ] )'] ;
                            prev.XDataSources{1} = ['meanNoNaN((',EpsApp,' - ',EpsF,') .* ', L0,'/', num2str(LcompMm{1}), '+ ',EpsF,',1)'] ;
                            prev.YDataSources{1} = 'hd.InputData-hd.InputData(1)' ;
                            prev.XDataSources{2} = ['meanNoNaN( hd.Seeds(',num2str(lstFil),').Strains(:,1,:,',num2str(ncam{1}),'),1)'] ;
                            prev.YDataSources{2} = 'hd.InputData-hd.InputData(1)' ;
                            prev.lines(1) = plot(NaN,NaN,'tag','Comp') ;
                            prev.lines(2) = plot(NaN,NaN,'tag','Fil') ;
                            prev.Axes.ColorOrderIndex = prev.Axes.ColorOrderIndex-1 ;
                            prev.legend = legend(prev.Axes,{'Comp','Fil'}) ;
                            prev.timeMarkers(1) = plot(NaN,NaN) ;
                            prev.timeMarkers(2) = plot(NaN,NaN) ;
                    end
                    title(regexprep(plotMachin,{'_'},{' vs '}))
                    prev = updatePreview(prev,hd) ;
                    col = colormap('lines') ;
                    for i = 1 : length(prev.lines)
                        set(prev.lines(i),'linestyle','-.','linewidth',1,'marker','.','markersize',20,'color',col(i,:))
                        set(prev.timeMarkers(i),'linestyle','none','marker','o','markersize',12,'linewidth',2,'color',col(i,:))
                        %if strcmpi(plotMachin,'FORCE_STRAIN3D')
                         %   set(prev.lines((i-1)*2+1),'linestyle','-.','linewidth',1,'marker','.','markersize',20,'color',col((i-1)*2+1,:))
                          %  set(prev.lines(i*2),'linestyle','-.','linewidth',1,'marker','.','markersize',20,'color',col(i*2,:))
                           % set(prev.timeMarkers((i-1)*2+1),'linestyle','none','marker','o','markersize',12,'linewidth',2,'color',col((i-1)*2+1,:))
                            %set(prev.timeMarkers(i*2),'linestyle','none','marker','o','markersize',12,'linewidth',2,'color',col(i*2,:))
                        %end
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
                                xdata = squeeze(eval(prev.XDataSources{l})) ;
                                ydata = squeeze(eval(prev.YDataSources{l})) ;
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
            
            function prev = listeSource(prev,hd)
                % Superclass updating function
                    prev.sources = {hd.Seeds(:).Name} ;
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