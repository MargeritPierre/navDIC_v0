function startInputAcquisitionFunction()

    answer = questdlg('Are you sure you want to take new datas ?', 'Connexion','Yes','No','Yes') ;
    if ~ismember(answer,'Yes') ; return ; end
% TIMER FOR PLOTTING
%     global plotVar;
%     plotVar = [];
%     plotVar.plotTimer = timer('ExecutionMode','FixedRate'...
%                     ,'Period',5 ...
%                     ,'TimerFcn',@(src,evt)plotTimerFunction());

% SIGNALS ACQUISITION
    global hd
%     
%     if exist('session','var') && isvalid(session)
%         session.NotifyWhenDataAvailableExceeds = 1 ; %to flush all undelivered datas
%         stop(session);
%     else
%         hd.DAQInputs.Session;
%     end

%     if exist('lst','var') && isvalid(lst) ; delete(lst) ; end
%     lst = addlistener(hd.DAQInputs.Session,'DataAvailable',@onDataAvailable);
    
    hd.DAQInputs.DataAcquisition.ScansAvailableFcn = @onDataAvailable ;
 
    hd.AcquiredData = [] ;
    hd.AcquiredData.Data = [] ;
    hd.AcquiredData.Time = [];
    hd.AcquiredData.StartTime = [];

    hd.DAQInputs.DataAcquisition.Rate = 800 ;
    %hd.DAQInputs.Session.TriggersPerRun = 1 ;
    %hd.DAQInputs.Session.IsContinuous = true ;
    hd.DAQInputs.DataAcquisition.ScansAvailableFcnCount = 8000 ;
    hd.AcquiredData.StartTime = now;
    
    start(hd.DAQInputs.DataAcquisition, 'Continuous');
    %startBackground(hd.DAQInputs.Session) ;
    
% MARKERS POSITION ACQUISITION

%     CreatePlot;  
    optTrack = hd.OT;
    
    hd.AcquiredOTPositions = [];
    hd.AcquiredOTPositions.Data = [];
    hd.AcquiredOTPositions.ID = [];
    hd.AcquiredOTPositions.Time = [];

    %if exist ('lst2','var') && isvalid(lst2) ; delete (lst2) ; end
    optTrack.addlistener(1, 'Position');
    hd.AcquiredOTPositions.StartTime = now;
    optTrack.enable(1);
%     
%     pos (= scatter(0, 0 ,500, [0 0.79 0.61] , '.');
    
%     start(plotVar.plotTimer);
    
    % CALLED FUNCTIONS
    function onDataAvailable(src,evt)
        [eventData, eventTimestamps, triggerTime] = read(src, src.ScansAvailableFcnCount, ...
            'OutputFormat', 'Matrix');
        hd.AcquiredData.Data = [hd.AcquiredData.Data ; eventData];
        disp(['Number of sample acquired: ', num2str(size(hd.AcquiredData.Data,1))])
        hd.AcquiredData.Time = [hd.AcquiredData.Time ; eventTimestamps] ;
        %evt.Source.read.Timestamps
    end

    % DISPLAY MARKERS POSITIONS
    function CreatePlot
        
        % create a figure which will contain two subplots
        f3 = figure('Name','Position of Optitrack Markers','NumberTitle','off');
        set(gca,'DataAspectRatio',[1,1,1]);
        %hf1.WindowStyle = 'docked';

        % plot and animated line for position
        %a1 = subplot( 1,2,1 );
        title( 'Markers Positions' );
        xlabel( 'x (m)' )
        ylabel( 'z (m)' )  
        
    end

% FUNCTION EXECUTED WHEN PLOT TIMER STARTS
    function plotTimerFunction( src, evnt)
        X = hd.AcquiredOTPositions.Data(:,3,end);
        Z = hd.AcquiredOTPositions.Data(:,2,end);
        set(pos,'XData', X, 'YData', Z)
        drawnow
    end

    
end
