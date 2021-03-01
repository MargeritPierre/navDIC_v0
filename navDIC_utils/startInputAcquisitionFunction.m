function startInputAcquisitionFunction()

    answer = questdlg('Are you sure you want to take new datas ?', 'Connexion','Yes','No','Yes') ;
    if ~ismember(answer,'Yes') ; return ; end
% TIMER FOR PLOTTING
    global plotVar;
    plotVar = [];
    plotVar.plotTimer = timer('ExecutionMode','FixedRate'...
                    ,'Period',1 ...
                    ,'TimerFcn',@(src,evt)plotTimerFunction());

% SIGNALS ACQUISITION
    global hd
    session = hd.DAQInputs.Session;

    if exist('lst','var') && isvalid(lst) ; delete(lst) ; end
    lst = addlistener(session,'DataAvailable',@onDataAvailable); 
 
    hd.AcquiredData = [] ;
    hd.AcquiredData.Data = [] ;
    hd.AcquiredData.Time = [];
    hd.AcquiredData.StartTime = [];

    session.Rate = 1000 ;
    session.TriggersPerRun = 1 ;
    session.IsContinuous = true ;
    session.NotifyWhenDataAvailableExceeds = 100 ;
    hd.AcquiredData.StartTime = now;
    startBackground(session) ;
    
% MARKERS POSITION ACQUISITION

    %CreatePlot;  
    optTrack = hd.OT;
    
    hd.AcquiredOTPositions = [];
    hd.AcquiredOTPositions.Data = [];
    hd.AcquiredOTPositions.ID = [];
    hd.AcquiredOTPositions.Time = [];

    %if exist ('lst2','var') && isvalid(lst2) ; delete (lst2) ; end
    optTrack.addlistener(1, 'Position');
    optTrack.enable(1);
    
    pos = scatter(0, 0 ,500, [0 0.79 0.61] , '.');
    
    start(plotVar.plotTimer);
    
    % CALLED FUNCTIONS
    function onDataAvailable(src,evt)
        hd.AcquiredData.Data = [hd.AcquiredData.Data ; evt.Data] ;
        hd.AcquiredData.Time = [hd.AcquiredData.Time ; evt.TimeStamps] ;
        %disp(['Number of sample acquired: ', num2str(size(hd.AcquiredData.Data,1))])
    end

    % DISPLAY MARKERS POSITIONS
    function CreatePlot
        
        % create a figure which will contain two subplots
        f1 = figure('Name','Position of Optitrack Markers','NumberTitle','off');
        set(gca,'DataAspectRatio',[1,1,1]);
        %hf1.WindowStyle = 'docked';

        % plot and animated line for position
        %a1 = subplot( 1,2,1 );
        title( 'Markers Positions' );
        xlabel( 'x (m)' )
        ylabel( 'z (m)' )  
        
    end

% FUNCTION EXECUTED WHEN PLOT TIMER STARTS
    function plotTimerFunction()
        X = hd.AcquiredOTPositions.Data(:,1,end);
        Z = hd.AcquiredOTPositions.Data(:,3,end);
        set(pos,'XData', X, 'YData', Z)
        drawnow
    end

    
end
