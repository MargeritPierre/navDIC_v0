function [Macros,macroListChanged] = manageMacros(Macros)
        
% ===================================================================================================================    
% MAIN FUNCTION
% ===================================================================================================================
    
    clc
    disp([newline,'MACRO MANAGEMENT :']) ;
    
    % Init Macro List
        availableMacros = [] ;
        getAvailableMacros() ;
        if isempty(availableMacros) ; return ; end
        if isempty(Macros) 
            newMacros = [] ; % Declare an empty structure
        else
            newMacros = copy(Macros) ;
        end
        
    % Init the figure
        fig = [] ;
        listBoxAvailable = [] ;
        listBoxCurrent = [] ;
        initFigure() ;
        updateLists() ;
        
    % Wait for the figure to be closed 
        macroListChanged = false ;
        validOutput = false ; % if not true, the function does not change anything (theoretically...)
        uiwait(fig) ;
        
    % Return Macros
        if validOutput % Changes will be made
            Macros = newMacros ;
        else % No changes
            macroListChanged = false ;
        end
        
        
% ===================================================================================================================    
% UTIL FUNCTIONS
% ===================================================================================================================

% UPDATE MACRO LISTS
    function updateLists()
        % ListBoxes on Figure
            listBoxAvailable.String = availableMacros ;
            if ~isempty(newMacros) 
                listBoxCurrent.String = {newMacros.Name} ; 
            else
                listBoxCurrent.String = {} ; 
            end
            listBoxAvailable.Value = 1 ;
            listBoxCurrent.Value = 1 ;
    end
            
% ADD MACRO
    function addMacro()
        % Get the selected item
            id = listBoxAvailable.Value ;
            macroToAdd = availableMacros{id} ;
        % Display the macro
            disp(newline) ;
            disp('NEW MACRO ADDED :') ;
            disp(macroToAdd) ;
        % Add the Macro
            macroToAdd = eval(macroToAdd) ;
            if isempty(newMacros) 
                newMacros = macroToAdd ;
            else
                newMacros(end+1) = macroToAdd ;
            end
        % Update Lists
            updateLists() ;
            macroListChanged = true ;
    end


% REMOVE MACRO
    function removeMacro()
        % If no used macros, return
            if isempty(newMacros) ; return ; end
        % Get the selected item
            id = listBoxCurrent.Value ;
            macroToRemove = newMacros(id) ;
        % Display a warning dialog
            answer = questdlg(['ARE YOU SURE YOU WANT TO REMOVE MACRO ' macroToRemove.Name '?']) ;
            if ~strcmpi(answer,'YES') ; return ; end
        % Display a message
            disp(newline) ;
            disp(['(' macroToRemove.Name ') MACRO HAS BEEN REMOVED']) ;
        % Remove the macro
            newMacros = newMacros(setdiff(1:length(newMacros),id)) ;
            delete(macroToRemove) ;
        % Update Lists
            updateLists() ;
            macroListChanged = true ;
    end

% MODIFY MACRO INFOS
    function modifyMacro()
        % If no used macros, return
            if isempty(newMacros) ; return ; end
        % Get the selected item
            id = listBoxCurrent.Value ;
        % Modify the macro settings
            global hd
            hd = setupUI(newMacros(id),hd) ;
        % Update Lists
            updateLists() ;
            macroListChanged = true ;
    end


% GET THE LIST OF AVAILABLE MACROS
    function getAvailableMacros()
        % Initialize
            availableMacros = [] ;
        % Get macro files
            global hd
            mfiles = dir([hd.RootPath,'/navDIC_macros/macros/*.m']) ;
            availableMacros = regexprep({mfiles.name},'\.m','') ;
    end


% OK BUTTON
    function clickOK()
        % Display a warning dialog
            if macroListChanged
                answer = questdlg(['ARE YOU SURE YOU WANT TO MAKE THE MODIFICATIONS ?']) ;
                if ~strcmpi(answer,'YES') ; return ; end
            end
        % Set the boolean to GO
            validOutput = true ;
        % Close the figure
            close(fig) ;
    end

% CANCEL BUTTON
    function clickCancel()
        % Display a warning dialog
            if macroListChanged
                answer = questdlg(['ARE YOU SURE YOU WANT TO UNDO THE MODIFICATIONS ?']) ;
                if ~strcmpi(answer,'YES') ; return ; end
            end
        % Set the boolean to GO
            validOutput = false ;
        % Close the figure
            close(fig) ;
    end

% CLEAR BUTTON
    function clickClear()
        % Display a warning dialog
            answer = questdlg(['ARE YOU SURE YOU WANT TO CLEAR THE MACROS ?']) ;
            if ~strcmpi(answer,'YES') ; return ; end
        % Update Lists
            updateLists() ;
            macroListChanged = true ;
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
            fig = figure('tag','manageMacros',...
                            'NumberTitle','off',...
                            'Name','MANAGE MACROS',...
                            'windowstyle','modal',...
                            'toolbar','none',...
                            'menubar','none') ;
            fig.Units = 'normalized' ;
            fig.Position = [.5-figSize/2 .5-figSize/2 figSize figSize] ;
        % First list
            listBoxAvailable = uicontrol(fig,'style','listbox') ;
            listBoxAvailable.Units = 'normalized' ;
            listBoxAvailable.Position = [margin 2*margin+btnHeight .5-2*margin-btnHeight/2 1-4*margin-btnHeight-ttlHeight] ;
        % First list Title
            ttl1 = uicontrol(fig,'style','text','string','AVAILABLE MACROS','BackgroundColor','w','fontweight','bold','units','normalized') ;
            ttl1.Position = [(listBoxAvailable.Position(1)+sum(listBoxAvailable.Position([1,3])))/2-ttlWidth/2 1-margin-ttlHeight ttlWidth ttlHeight] ;
        % Second list
            listBoxCurrent = uicontrol(fig,'style','listbox') ;
            listBoxCurrent.Units = 'normalized' ;
            listBoxCurrent.Position = [.5+1*margin+btnHeight/2 2*margin+btnHeight .5-2*margin-btnHeight/2 1-4*margin-btnHeight-ttlHeight] ;
        % Second list Title
            ttl2 = uicontrol(fig,'style','text','string','USED MACROS','BackgroundColor','w','fontweight','bold','units','normalized') ;
            ttl2.Position = [(listBoxCurrent.Position(1)+sum(listBoxCurrent.Position([1,3])))/2-ttlWidth/2 1-margin-ttlHeight ttlWidth ttlHeight] ;
        % AddMacro Button
            btnAdd = uicontrol('style','pushbutton','units','normalized','string','+','fontweight','bold') ;
            btnAdd.Position = [.5-btnHeight/2 .5+btnHeight/2+margin/2 btnHeight btnHeight] ;
            btnAdd.Callback = @(src,evt)addMacro() ;
        % RemoveMacro Button
            btnRemove = uicontrol('style','pushbutton','units','normalized','string','-','fontweight','bold') ;
            btnRemove.Position = [.5-btnHeight/2 .5-btnHeight/2 btnHeight btnHeight] ;
            btnRemove.Callback = @(src,evt)removeMacro() ;
        % ModifyInfos Button
            btnInfos = uicontrol('style','pushbutton','units','normalized','string','Setup','fontweight','bold') ;
            btnInfos.Position = [.5-btnHeight/2 .5-3*btnHeight/2-margin btnHeight btnHeight] ;
            btnInfos.Callback = @(src,evt)modifyMacro() ;
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