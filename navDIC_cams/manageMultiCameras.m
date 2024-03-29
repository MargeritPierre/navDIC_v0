function [CAMERAS,camsHasChanged] = manageMultiCameras(CAMERAS)
        
% ===================================================================================================================    
% MAIN FUNCTION
% ===================================================================================================================
    
    clc
    display([char(10),'CAMERAS MANAGEMENT UTIL :']) ;
    
    % Init IMAQ Session
        if isempty(CAMERAS) 
            imaqreset ;
            CAMERAS = [] ;
        end
        camsHasChanged = false ;
    
    % Init CAMERA List
        adaptors = [] ;
        availableCams = [] ;
        getAvailableCams() ;
        if isempty(availableCams) ; return ; end
        usedCams = [] ;
        freeCams = [] ;
        % UsedInputs
            if isempty(CAMERAS) 
                usedCams = [] ; % Declare an empty structure
            else
                usedCams = CAMERAS ;
            end 
        
    % Init the figure
        fig = [] ;
        listBoxFree = [] ;
        listBoxUsed = [] ;
        initFigure() ;
        updateLists() ;
        
    % Wait for the figure to be closed 
        while ishandle(fig) ; drawnow ; end
        
    % Return CAMERAS
        CAMERAS = usedCams ;
        
        
% ===================================================================================================================    
% UTIL FUNCTIONS
% ===================================================================================================================

% SET HARDWARE INPUT
    function cam = connectCamera(cam)
        % Choose the Video Format
            DeviceInfos = imaqhwinfo(cam.Adaptor, cam.Infos.DeviceID) ;
            [id,valid] = listdlg('PromptString','Select a Video Format:',...
                                'SelectionMode','single',...
                                'initialValue',1,...
                                'ListString',DeviceInfos.SupportedFormats) ;
            if ~valid ; return ; end
        % Declare the videoinput
            cam.VidObj = videoinput(cam.Adaptor, cam.Infos.DeviceID,DeviceInfos.SupportedFormats{id}) ;
        % Reglage du trigger
            listTriggerType = {triggerinfo(cam.VidObj).TriggerType} ;
            %listTriggerType = [{'manual'}, {'hardware'}, {'infinite'}] ;
            [id,valid] = listdlg('PromptString','Select a trigger type :',...
                'SelectionMode','single',...
                'initialValue',1,...
                'ListString',listTriggerType) ;
            if ~valid ; return ; end
            if strcmp(cam.VidObj.Running,'on')
                stop(cam) ;   
            end
            triggerconfig(cam.VidObj,listTriggerType{id}) 
        % Set RBG ColorSpace
            if ~strcmpi(cam.VidObj.ReturnedColorSpace,'grayscale')
                cam.VidObj.ReturnedColorSpace = 'RGB' ;
            end
        % Set trigger mode            
            cam.VidObj.FramesPerTrigger = 1 ;
            cam.VidObj.TriggerRepeat = Inf ;
        % Set FrameAcquiredFunction
            cam.VidObj.FramesAcquiredFcnCount = 1 ;
            cam.VidObj.FramesAcquiredFcn = @(src,evt)flushdata(src) ;
        % Retrieve infos
            cam.Infos.IMAQ = propinfo(cam.VidObj.Source) ;
        % Add custom informations
            cam.CurrentState = 'connected' ;
    end

% RESET ALL HARDWARE INPUTS
    function resetSession()
        % Delete All Cameras
            delete(imaqfind) ;
            imaqreset ;
            delete(imaqfind) ;
        % Re-initialize Lists
            usedCams = [] ;
            getAvailableCams() ;
            updateLists() ;
    end

% UPDATE INPUT LISTS
    function updateLists()
        % Free Cameras
            if isempty(usedCams) 
                freeCams  = availableCams ;
            else
                isUsed = ismember({availableCams.Name},{usedCams.Name}) ;
                freeCams = availableCams(~isUsed) ;
            end
        % ListBoxes on Figure
            if ~isempty(freeCams)
                listBoxFree.String = {freeCams.Name} ;
            else
                listBoxFree.String = {} ; 
            end
            if ~isempty(usedCams) 
                listBoxUsed.String = {usedCams.Name} ; 
            else
                listBoxUsed.String = {} ; 
            end
            listBoxFree.Value = 1 ;
            listBoxUsed.Value = 1 ;
    end
            
% ADD INPUT
    function addCamera()
        % Get the selected item
            id = listBoxFree.Value ;
            camToAdd = freeCams(id) ;
        % Launch the Hardware Connection
            camToAdd = connectCamera(camToAdd) ;
            if ~isfield(camToAdd,'VidObj') ; return ; end
            %start(camToAdd.VidObj) ; %<NEW>
        % Display the input
            disp(char(10)) ;
            disp('NEW INPUT ADDED :') ;
            disp(camToAdd) ;
        % Add the Input
            if isempty(usedCams) 
                usedCams = camToAdd ;
            else
                usedCams(end+1) = camToAdd ;
            end
        % Update Lists
            updateLists() ;
            camsHasChanged = true ;
    end


% REMOVE INPUT
    function removeCamera()
        % If no used inputs, return
            if isempty(usedCams) ; return ; end
        % Get the selected item
            id = listBoxUsed.Value ;
            camToRemove = usedCams(id) ;
        % Display a warning dialog
            answer = questdlg(['ARE YOU SURE YOU WANT TO REMOVE CAMERA ' camToRemove.Infos.DeviceName '?'],'DELETE CAMERA ?','OK','Cancel','Cancel') ;
            if ~strcmp(upper(answer),'OK') ; return ; end
        % Release the camera
            delete(camToRemove.VidObj) ;
            delete(usedCams(id).VidObj) ;
        % Removed the item in the list
            usedCams = usedCams(setdiff(1:length(usedCams),id)) ;
        % Display a message
            disp(char(10)) ;
            disp([camToRemove.Infos.DeviceName ' CAMERA HAS BEEN REMOVED']) ;
        % Update Lists
            updateLists() ;
            camsHasChanged = true ;
    end

% MODIFY INPUT INFOS
    function modifyCameraSettings()
        % If no used inputs, return
            if isempty(usedCams) ; return ; end
        % Get the selected item
            id = listBoxUsed.Value ;
        % Re-Set Infos
            updatedCam = setCameraSettings(usedCams(id)) ;
            usedCams(id) = updatedCam ;
        % Update Lists
            updateLists() ;
            camsHasChanged = true ;
    end


% GET THE LIST OF AVAILABLE INPUTS
    function getAvailableCams()
        % Initialize
            adaptors = [] ;
            availableCams = [] ;
            infos = imaqhwinfo ;
        % Get installed adaptors
            adaptors = infos.InstalledAdaptors ;
        % Connected cameras
            nA = length(adaptors) ;
            if ~nA ; disp('NO INSTALLED ADAPTORS') ; return ; end
            for a = 1:nA
                adaptName = adaptors{a} ;
                disp(['   ',adaptName]) ;
                % Cameras connected with this adaptor
                    adaptCams = imaqhwinfo(adaptName) ;
                    nC = length(adaptCams.DeviceIDs) ;
                    if ~nC ; disp('      NO CONNECTED CAMERAS') ; continue ; end
                        for c = 1:nC
                            availableCams(end+1).Infos = adaptCams.DeviceInfo(c) ;
                            availableCams(end).Adaptor = adaptName ;
                            availableCams(end).Name = [availableCams(end).Infos.DeviceName ' ' num2str(availableCams(end).Infos.DeviceID)] ;
                        end
            end
    end


% OK BUTTON
    function clickOK()
        % Close the figure
            fig.CloseRequestFcn = @(src,evt)closereq ;
            close(fig) ;
    end

% CANCEL BUTTON
    function clickCancel()
    end

% CLEAR BUTTON
    function clickClear()
        % Is there some Cams to clear ?
            if ~isempty(usedCams)
                % Display a warning dialog
                    answer = questdlg(['ARE YOU SURE YOU WANT TO CLEAR ALL CAMERAS ?'],'CLEAR ALL CAMERAS ?','OK','Cancel','Cancel') ;
                    if ~strcmp(upper(answer),'OK') ; return ; end
                % Reset the session
                    resetSession() ;
                    camsHasChanged = true ;
            end
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
                            'Name','MANAGE CAMERAS',...
                            'windowstyle','modal',...
                            'toolbar','none',...
                            'menubar','none',...
                            'NumberTitle','off'...
                            ) ;
            fig.Units = 'normalized' ;
            fig.Position = [.5-figSize/2 .5-figSize/2 figSize figSize] ;
            fig.CloseRequestFcn = @(src,evt)clickOK ;
        % First list
            listBoxFree = uicontrol(fig,'style','listbox') ;
            listBoxFree.Units = 'normalized' ;
            listBoxFree.Position = [margin 2*margin+btnHeight .5-2*margin-btnHeight/2 1-4*margin-btnHeight-ttlHeight] ;
        % First list Title
            ttl1 = uicontrol(fig,'style','text','string','FREE CAMERAS','BackgroundColor','w','fontweight','bold','units','normalized') ;
            ttl1.Position = [(listBoxFree.Position(1)+sum(listBoxFree.Position([1,3])))/2-ttlWidth/2 1-margin-ttlHeight ttlWidth ttlHeight] ;
        % Second list
            listBoxUsed = uicontrol(fig,'style','listbox') ;
            listBoxUsed.Units = 'normalized' ;
            listBoxUsed.Position = [.5+1*margin+btnHeight/2 2*margin+btnHeight .5-2*margin-btnHeight/2 1-4*margin-btnHeight-ttlHeight] ;
        % Second list Title
            ttl2 = uicontrol(fig,'style','text','string','USED CAMERAS','BackgroundColor','w','fontweight','bold','units','normalized') ;
            ttl2.Position = [(listBoxUsed.Position(1)+sum(listBoxUsed.Position([1,3])))/2-ttlWidth/2 1-margin-ttlHeight ttlWidth ttlHeight] ;
        % AddInput Button
            btnAdd = uicontrol('style','pushbutton','units','normalized','string','+','fontweight','bold') ;
            btnAdd.Position = [.5-btnHeight/2 .5+btnHeight/2+margin/2 btnHeight btnHeight] ;
            btnAdd.Callback = @(src,evt)addCamera() ;
        % RemoveInput Button
            btnRemove = uicontrol('style','pushbutton','units','normalized','string','-','fontweight','bold') ;
            btnRemove.Position = [.5-btnHeight/2 .5-btnHeight/2 btnHeight btnHeight] ;
            btnRemove.Callback = @(src,evt)removeCamera() ;
        % ModifyInfos Button
            btnInfos = uicontrol('style','pushbutton','units','normalized','string','Settings','fontweight','bold') ;
            btnInfos.Position = [.5-btnHeight/2 .5-3*btnHeight/2-margin btnHeight btnHeight] ;
            btnInfos.Callback = @(src,evt)modifyCameraSettings() ;
        % OK Button
            btnWidth = (1-4*margin)/3 ;
            btnOK = uicontrol('style','pushbutton','units','normalized','string','OK','fontweight','bold') ;
            btnOK.Position = [margin margin btnWidth btnHeight] ;
            btnOK.Callback = @(src,evt)clickOK() ;
        % CANCEL Button
            btnCancel = uicontrol('style','pushbutton','units','normalized','string','Cancel','fontweight','bold') ;
            btnCancel.Position = [2*margin+btnWidth margin btnWidth btnHeight] ;
            btnCancel.Callback = @(src,evt)clickCancel() ;
            btnCancel.Enable = 'off' ;
        % OK Button
            btnClear = uicontrol('style','pushbutton','units','normalized','string','Clear','fontweight','bold') ;
            btnClear.Position = [2*margin+2*btnWidth margin btnWidth btnHeight] ;
            btnClear.Callback = @(src,evt)clickClear() ;
    end





end