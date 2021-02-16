function [OTModel, OTinputsListChanged] = manageOTInputs()

% ===================================================================================================================    
% MAIN FUNCTION
% ===================================================================================================================

    % Create an instance of the natnet client class
        OT = natnet;

    % Init the figure
        fig = [] ;
        initFigure() ;
        
        % Wait for the figure to be closed 
        OTinputsListChanged = false ;      
        uiwait(fig) ;

    % Set Motive in Live Mode for record
    if (OT.Mode == 'Edit')
        OT.liveMode;
    end
    
% ===================================================================================================================    
% UTIL FUNCTIONS
% ===================================================================================================================    
  

% ConnectOT Button
    function ConnectOT()
        % Connect the client to the server (multicast over local loopback)
        fprintf('Connecting to the server ...\n')
        OT.HostIP = '127.0.0.1';
        OT.ClientIP = '127.0.0.1';
        OT.ConnectionType = 'Multicast';
        OT.connect;
        if ( OT.IsConnected == 0 )
            fprintf( 'Failed to connect\n' )
            fprintf( '\tMake sure the host is connected to the network\n' )
            fprintf( '\tand that the host and client IP addresses are correct\n\n' ) 
        return
        else
            fprintf( 'Successfully connected to OptiTracks\n' )
            OTModel = OT.getModelDescription;
            if isempty(OTModel)
                fprintf( 'Could not get the model. Did you calibrated the cameras ?\n' )
            end
        end
        % Close the figure
        close(fig); 
    end

% Cancel Button
    function clickCancel()
        close(fig) ;
    end
    
% Init the figure

    function initFigure()
        % Parameters (all sizes are normalized)
            figWidth = 0.25 ;
            figHeigth = .15 ;
            margin = .05 ;
            btnHeight = .15 ;
            btnWidth = .2 ;
            ttlHeight = .2 ;
            ttlWidth = .6 ;
        % Figure centered on the screen
            fig = figure('tag','manageOTInputs',...
                            'NumberTitle','off',...
                            'Name','MANAGE OPTITRACKS INPUTS',...
                            'windowstyle','modal',...
                            'toolbar','none',...
                            'menubar','none') ;
            fig.Units = 'normalized' ;
            fig.Position = [.5-figWidth/2 .5-figHeigth/2 figWidth figHeigth] ;
        % Text
            ttl1 = uicontrol('style','text','string','Do you want to connect to Motive ?','BackgroundColor','w','fontweight','bold','FontSize',10,'units','normalized') ;
            ttl1.Position = [.5-figWidth-margin .5+figHeigth ttlWidth ttlHeight] ;            
        % AddConnectOT Button
            btnOK = uicontrol('style','pushbutton','units','normalized','string','YES','fontweight','bold') ;
            btnOK.Position = [.5-1.2*btnWidth .5-figHeigth btnWidth btnHeight] ;
            btnOK.Callback = @(src,evt)ConnectOT() ;
        % CANCEL Button
            btnOK = uicontrol('style','pushbutton','units','normalized','string','NO','fontweight','bold') ;
            btnOK.Position = [.5+0.1*btnWidth .5-figHeigth btnWidth btnHeight] ;
            btnOK.Callback = @(src,evt)clickCancel() ;
    end

end