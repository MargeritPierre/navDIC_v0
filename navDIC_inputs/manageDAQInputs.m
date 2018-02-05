function [DAQInputs,inputListChanged] = manageDAQInputs(DAQInputs)
        
% ===================================================================================================================    
% MAIN FUNCTION
% ===================================================================================================================
    
    clc
    display([char(10),'INPUT MANAGEMENT :']) ;
    
    % Init DAQ Session
        if ~isfield(DAQInputs,'Session') 
            daqreset ;
            DAQInputs.Session = daq.createSession('ni') ;
            DAQInputs.Inputs = {} ;
        end
    
    % Init Input List
        devices = [] ;
        availableInputs = [] ;
        getAvailableInputs() ;
        if isempty(availableInputs) ; return ; end
        usedInputs = [] ;
        freeInputs = [] ;
        % UsedInputs
            if isempty(DAQInputs.Inputs) 
                usedInputs = [] ; % Declare an empty structure
            else
                usedInputs = DAQInputs.Inputs ;
            end
        
    % Init the figure
        fig = [] ;
        listBoxFree = [] ;
        listBoxUsed = [] ;
        initFigure() ;
        updateLists() ;
        
    % Wait for the figure to be closed 
        inputListChanged = false ;
        validOutput = false ; % if not true, the function does not change anything (theoretically...)
        uiwait(fig) ;
        
    % Return DAQInputs
        if validOutput % Changes will be made
            DAQInputs.Inputs = usedInputs ;
        else % No changes
            DAQInputs.Inputs = resetSession(DAQInputs.Inputs) ;
            inputListChanged = false ;
        end
        
        
% ===================================================================================================================    
% UTIL FUNCTIONS
% ===================================================================================================================

% SET HARDWARE INPUT
    function input = setHardInput(input)
        input.HardwareInput = ...
            addAnalogInputChannel(...
            DAQInputs.Session,...
            devices(input.DeviceID).ID,...
            devices(input.DeviceID).Subsystems(input.SubSystID).ChannelNames{input.ChannelID},...
            'Voltage') ;
        input.HardwareInput.Range = input.Ranges(input.RangeID) ;
    end

% SET INPUT INFO
    function [input,valid] = setInputInfos(input)
        % Get AvailableRanges
            ranges = devices(input.DeviceID).Subsystems(input.SubSystID).RangesAvailable ;
            nR = length(ranges) ;
            rangeNames = {} ;
            indexedRangeNames = {} ;
            for r = 1:nR
                rangeNames{end+1} = [num2str(ranges(r).Min) ' to ' num2str(ranges(r).Max) ' ' ranges(r).Units] ;
                indexedRangeNames{end+1} = [num2str(r) ' : ',rangeNames{end}] ;
            end
            joinedRangeNames = strjoin(indexedRangeNames,char(10)) ;
        % Get User-Defined Input Infos
            valid = false ;
            askFor = {'NAME','UNITS','SENSITIVITY (Units/Volts)',['RANGE (TYPE THE ID)' char(10) joinedRangeNames]} ;
            title = 'USER-DEFINED INPUT INFOS';
            num_lines = 1 ;
            % Default Values
                if isfield(input,'DataName') % Not the first time infos are setted
                    default = {input.DataName,...
                                input.Units,...
                                num2str(input.Sensitivity),...
                                num2str(input.RangeID)};
                else
                    default = {['in',num2str(length(usedInputs)+1)],'V','1','1'};
                end
            answers = inputdlg(askFor,title,num_lines,default) ;
            if isempty(answers) ; return ; end
            % Extract answers
                input.DataName = answers{1} ;
                input.Units = answers{2} ;
                input.Sensitivity = str2num(answers{3}) ;
                input.RangeID = str2num(answers{4}) ;
                input.Ranges = ranges ;
        % Set the FullName
            input.FullName = [...
                                    input.DataName ...
                                    ', ' ...
                                    num2str(input.Sensitivity) ...
                                    ' ' ...
                                    input.Units ...
                                    '/V ' ...
                                    ', ' ...
                                    rangeNames{input.RangeID} ...
                                    ' (' ...
                                    input.Name ...
                                    ')'...
                                    ] ;
        % The input infos are valid
            valid = true ;
    end

% RESET ALL HARDWARE INPUTS
    function keptInputs = resetSession(keptInputs)
        % Release the Session
            release(DAQInputs.Session) ;
        % Re-create a Session
            DAQInputs.Session = daq.createSession('ni') ;
        % And re-set ALL harware inputs
            if ~isempty(keptInputs)
                for in = 1:length(keptInputs)
                    keptInputs(in) = setHardInput(keptInputs(in)) ;
                end
            end
    end

% UPDATE INPUT LISTS
    function updateLists()
        % Free Inputs
            if isempty(usedInputs) 
                freeInputs  = availableInputs ;
            else
                [~,freeInputsIndices] = setdiff({availableInputs.Name},{usedInputs.Name}) ;
                freeInputs = availableInputs(freeInputsIndices) ;
            end
        % ListBoxes on Figure
            if ~isempty(freeInputs)
                listBoxFree.String = {freeInputs.Name} ;
            else
                listBoxFree.String = {} ; 
            end
            if ~isempty(usedInputs) 
                listBoxUsed.String = {usedInputs.FullName} ; 
            else
                listBoxUsed.String = {} ; 
            end
            listBoxFree.Value = 1 ;
            listBoxUsed.Value = 1 ;
    end
            
% ADD INPUT
    function addInput()
        % Get the selected item
            id = listBoxFree.Value ;
            inputToAdd = freeInputs(id) ;
        % Set Input Infos
            [inputToAdd,valid] = setInputInfos(inputToAdd) ;
            if ~valid ; return ; end
        % Set the Hardware Input
            inputToAdd = setHardInput(inputToAdd) ;
        % Display the input
            disp(char(10)) ;
            disp('NEW INPUT ADDED :') ;
            disp(inputToAdd) ;
        % Add the Input
            if isempty(usedInputs) 
                usedInputs = inputToAdd ;
            else
                usedInputs(end+1) = inputToAdd ;
            end
        % Update Lists
            updateLists() ;
            inputListChanged = true ;
    end


% REMOVE INPUT
    function removeInput()
        % If no used inputs, return
            if isempty(usedInputs) ; return ; end
        % Get the selected item
            id = listBoxUsed.Value ;
            inputToRemove = usedInputs(id) ;
        % Display a warning dialog
            answer = questdlg(['ARE YOU SURE YOU WANT TO REMOVE INPUT ' inputToRemove.DataName '?']) ;
            if ~strcmp(upper(answer),'YES') ; return ; end
        % One cannot release a single input... 
        % Removed the input
            usedInputs = usedInputs(setdiff(1:length(usedInputs),id)) ;
        % Release the ENTIRE Session
            usedInputs = resetSession(usedInputs) ;
        % Display a message
            disp(char(10)) ;
            disp(['(' inputToRemove.DataName ') INPUT HAS BEEN REMOVED']) ;
        % Update Lists
            updateLists() ;
            inputListChanged = true ;
    end

% MODIFY INPUT INFOS
    function modifyInputInfos()
        % If no used inputs, return
            if isempty(usedInputs) ; return ; end
        % Get the selected item
            id = listBoxUsed.Value ;
        % Re-Set Infos
            [updatedInput,valid] = setInputInfos(usedInputs(id)) ;
            if ~valid ; return ; end
            usedInputs(id) = updatedInput ;
        % Update Lists
            updateLists() ;
            inputListChanged = true ;
    end


% GET THE LIST OF AVAILABLE INPUTS
    function getAvailableInputs()
        % Initialize
            devices = [] ;
            availableInputs = [] ;
        % Get available devices
            devices = daq.getDevices() ;
        % Available devices
            nD = length(devices);
            if ~nD ; disp('NO AVAILABLE DEVICES') ; return ; end
            for d = 1:nD
                %devName = [dev(d).Vendor.FullName,', ',dev(d).Model] ;
                devName = devices(d).Model ;
                disp(['   ',devName]) ;
                % Subsystems of the device
                    nS = length(devices(d).Subsystems) ;
                    if ~nS ; disp('      NO AVAILABLE SUBSYSTEMS') ; return ; end
                    for s = 1:nS
                        subSystName = devices(d).Subsystems(s).SubsystemType ;
                        disp(['      ',subSystName])
                        % Inputs Available
                            nC = length(devices(d).Subsystems(s).ChannelNames) ;
                            if ~nC ; disp('         NO AVAILABLE INPUTS') ; return ; end
                            for c = 1:nC
                                channelName = devices(d).Subsystems(s).ChannelNames{c} ;
                                disp(['         ',channelName])
                                % Is this Channel an input ?
                                    if regexp(channelName,'ai')
                                        availableInputs(end+1).Name = [devName,', ',channelName] ;
                                        availableInputs(end).DeviceID = d ;
                                        availableInputs(end).SubSystID = s ;
                                        availableInputs(end).ChannelID = c ;
                                    end
                            end
                    end
            end
    end


% OK BUTTON
    function clickOK()
        % Display a warning dialog
            if inputListChanged
                answer = questdlg(['ARE YOU SURE YOU WANT TO MAKE THE MODIFICATIONS ?']) ;
                if ~strcmp(upper(answer),'YES') ; return ; end
            end
        % Set the boolean to GO
            validOutput = true ;
        % Close the figure
            close(fig) ;
    end

% CANCEL BUTTON
    function clickCancel()
        % Display a warning dialog
            if inputListChanged
                answer = questdlg(['ARE YOU SURE YOU WANT TO UNDO THE MODIFICATIONS ?']) ;
                if ~strcmp(upper(answer),'YES') ; return ; end
            end
        % Set the boolean to GO
            validOutput = false ;
        % Close the figure
            close(fig) ;
    end

% CLEAR BUTTON
    function clickClear()
        % Display a warning dialog
            answer = questdlg(['ARE YOU SURE YOU WANT TO CLEAR THE SESSION ?']) ;
            if ~strcmp(upper(answer),'YES') ; return ; end
        % Reset the session with no keptInputs
            usedInputs = resetSession([]) ;
        % Update Lists
            updateLists() ;
            inputListChanged = true ;
    end


% INIT THE FIGURE
    function initFigure()
        % Parameters (all sizes are normalized)
            figSize = .3 ;
            margin = .01 ;
            btnHeight = .1 ;
            ttlHeight = .05 ;
            ttlWidth = .45 ;
        % Figure centered on the screen
            fig = figure('tag','manageDAQInputs',...
                            'NumberTitle','off',...
                            'Name','MANAGE DAQ INPUTS',...
                            'windowstyle','modal',...
                            'toolbar','none',...
                            'menubar','none') ;
            fig.Units = 'normalized' ;
            fig.Position = [.5-figSize/2 .5-figSize/2 figSize figSize] ;
        % First list
            listBoxFree = uicontrol(fig,'style','listbox') ;
            listBoxFree.Units = 'normalized' ;
            listBoxFree.Position = [margin 2*margin+btnHeight .5-2*margin-btnHeight/2 1-4*margin-btnHeight-ttlHeight] ;
        % First list Title
            ttl1 = uicontrol(fig,'style','text','string','FREE INPUTS','BackgroundColor','w','fontweight','bold','units','normalized') ;
            ttl1.Position = [(listBoxFree.Position(1)+sum(listBoxFree.Position([1,3])))/2-ttlWidth/2 1-margin-ttlHeight ttlWidth ttlHeight] ;
        % Second list
            listBoxUsed = uicontrol(fig,'style','listbox') ;
            listBoxUsed.Units = 'normalized' ;
            listBoxUsed.Position = [.5+1*margin+btnHeight/2 2*margin+btnHeight .5-2*margin-btnHeight/2 1-4*margin-btnHeight-ttlHeight] ;
        % Second list Title
            ttl2 = uicontrol(fig,'style','text','string','USED INPUTS','BackgroundColor','w','fontweight','bold','units','normalized') ;
            ttl2.Position = [(listBoxUsed.Position(1)+sum(listBoxUsed.Position([1,3])))/2-ttlWidth/2 1-margin-ttlHeight ttlWidth ttlHeight] ;
        % AddInput Button
            btnAdd = uicontrol('style','pushbutton','units','normalized','string','+','fontweight','bold') ;
            btnAdd.Position = [.5-btnHeight/2 .5+btnHeight/2+margin/2 btnHeight btnHeight] ;
            btnAdd.Callback = @(src,evt)addInput() ;
        % RemoveInput Button
            btnRemove = uicontrol('style','pushbutton','units','normalized','string','-','fontweight','bold') ;
            btnRemove.Position = [.5-btnHeight/2 .5-btnHeight/2 btnHeight btnHeight] ;
            btnRemove.Callback = @(src,evt)removeInput() ;
        % ModifyInfos Button
            btnInfos = uicontrol('style','pushbutton','units','normalized','string','Infos','fontweight','bold') ;
            btnInfos.Position = [.5-btnHeight/2 .5-3*btnHeight/2-margin btnHeight btnHeight] ;
            btnInfos.Callback = @(src,evt)modifyInputInfos() ;
        % OK Button
            btnWidth = (1-4*margin)/3 ;
            btnOK = uicontrol('style','pushbutton','units','normalized','string','OK','fontweight','bold') ;
            btnOK.Position = [margin margin btnWidth btnHeight] ;
            btnOK.Callback = @(src,evt)clickOK() ;
        % CANCEL Button
            btnOK = uicontrol('style','pushbutton','units','normalized','string','Cancel','fontweight','bold') ;
            btnOK.Position = [2*margin+btnWidth margin btnWidth btnHeight] ;
            btnOK.Callback = @(src,evt)clickCancel() ;
        % OK Button
            btnOK = uicontrol('style','pushbutton','units','normalized','string','Clear','fontweight','bold') ;
            btnOK.Position = [2*margin+2*btnWidth margin btnWidth btnHeight] ;
            btnOK.Callback = @(src,evt)clickClear() ;
    end





end