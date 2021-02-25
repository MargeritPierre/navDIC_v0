function [OT, nMarkers] = manageOTInputs()

% ===================================================================================================================    
% MAIN FUNCTION
% ===================================================================================================================
      
    % Create an instance of the natnet client class
        OT = natnet;

    % Init the figure
        fig = [] ;
        answer = questdlg('Do you want to connect to Motive ?', 'Connexion','Yes','No','Yes') ;
        if ~ismember(answer,'Yes') ; return ; end
        ConnectOT()
        
        nMarkers = 0 ; 
        if ~isempty(OT)
            for i = 1:1000
                if ~isempty(OT.getFrame.LabeledMarker(i))
                    nMarkers = nMarkers + 1 ;
                end
            end
        end
        %hd.AcquiredOTPositions.MarkersCount = nMarkers;
        
        % Wait for the figure to be closed 
        %OTinputsListChanged = false ;      
       %uiwait(fig) ;

    % Set Motive in Live Mode for record
    if ~strcmp(OT.Mode,'Live')
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

end