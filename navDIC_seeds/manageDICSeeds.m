function hd = manageDICSeeds(hd)

% ===================================================================================================================    
% MAIN FUNCTION
% ===================================================================================================================
    
    clc
    display([char(10),'SEED MANAGEMENT :']) ;
    
    % Init Seeds
        if ~isfield(hd,'Seeds') 
            hd.Seeds = [] ;
        end
    
    % Init Seed Classes List
        availableSeedClasses = [] ;
        getAvailableSeedClasses() ;
        if isempty(availableSeedClasses) 
            disp('NO SEED CLASSES AVAILABLE !') ;
            return ; 
        end
        SEEDS = hd.Seeds ;
        
    % Init the figure
        fig = [] ;
        popupClasses = [] ;
        listBoxSeeds = [] ;
        infoText = [] ;
        initFigure() ;
        updateSeedList() ;
        
    % Wait for the figure to be closed 
        inputListChanged = false ;
        validOutput = false ; % if not true, the function does not change anything (theoretically...)
        while isvalid(fig)
            drawnow ;
        end
        
    % Return DAQInputs
        if validOutput % Changes will be made
            hd.Seeds = SEEDS ;
        end
        
        
% ===================================================================================================================    
% UTIL FUNCTIONS
% ===================================================================================================================


% CREATE SEED
    function [seed,valid] = createSeed(seedClass)
        % Input a name for the seed
            newSeedName = inputdlg('NAME OF THE SEED ?','Input a Name',1,{seedClass.ClassName}) ;
            if isempty(newSeedName) ; return ; end
        % Initialize the new seed
            split = strsplit(seedClass.ClassName,'_') ;
            if strcmpi(split{1},'3D')
                if ~isfield(hd.Cameras,'Properties')
                    for i = 1 : length(hd.Cameras)
                        hd = setCameraProperties(hd,i) ;
                    end
                end
            end
            seed = eval([seedClass.Class,'(hd)']) ; %navDICSeed_2D_DistMesh(hd) ; %
            seed.Name = newSeedName{1} ;
            
        % Output OK
            valid = seed.isValid ;
    end

% MODIFY SEED
    function modifySeed()
        ID = listBoxSeeds.Value ;
        seed = SEEDS(ID) ;
        seed = seed.modify(hd) ;
        SEEDS(ID) = seed ;
        % Update Lists
            updateSeedList() ;
            inputListChanged = true ;
    end

% RESET ALL SEEDS
    function resetSeeds()
        SEEDS = [] ;
        getAvailableSeedClasses() ;
    end

% UPDATE INPUT LISTS
    function updateSeedList()
        % ListBox on Figure
            if ~isempty(SEEDS) 
                listBoxSeeds.String = {SEEDS.Name} ; 
            else
                listBoxSeeds.String = {} ; 
            end
            listBoxSeeds.Value = 1 ;
    end
            
% ADD SEED
    function addSeed()
        % Get the selected item
            id = popupClasses.Value ;
            seedClassToAdd = availableSeedClasses(popupClasses.Value) ;
        % Create the Seed
            [seedToAdd,valid] = createSeed(seedClassToAdd) ;
            if ~valid ; return ; end
        % Give a name to the seed
        % Display the input
            disp(char(10)) ;
            disp('NEW SEED ADDED :') ;
            disp(seedToAdd) ;
        % Add the Seed
            if isempty(SEEDS) 
                SEEDS = seedToAdd ;
            else
                SEEDS(end+1) = seedToAdd ;
            end
        % Update Lists
            updateSeedList() ;
            inputListChanged = true ;
    end


% REMOVE INPUT
    function removeSeed()
        % If no defined seeds, return
            if isempty(SEEDS) ; return ; end
        % Get the selected item
            id = listBoxSeeds.Value ;
            seedToRemove = SEEDS(id) ;
        % Display a warning dialog
            answer = questdlg(['ARE YOU SURE YOU WANT TO REMOVE SEED ' seedToRemove.Name '?']) ;
            if ~strcmp(upper(answer),'YES') ; return ; end
        % Remove the seed
            SEEDS = SEEDS(setdiff(1:length(SEEDS),id)) ;
        % Display a message
            disp(char(10)) ;
            disp(['(' seedToRemove.Name ') INPUT HAS BEEN REMOVED']) ;
        % Update Lists
            updateSeedList() ;
            inputListChanged = true ;
    end

% DUPLICATE A SEED
    function duplicateSeed()
        % If no defined seeds, return
            if isempty(SEEDS) ; return ; end
        % Get the selected item
            id = listBoxSeeds.Value ;
        % Input a name for the seed
            newSeedName = inputdlg('NAME OF THE SEED COPY ?','Input a Name',1,{['Copy_of_',SEEDS(id).Name]}) ;
            if isempty(newSeedName) ; return ; end
        % Add the seed copy
            SEEDS(end+1) = SEEDS(id) ;
            SEEDS(end).Name = newSeedName{1} ;
        % Update infos
            updateSeedList() ;
            inputListChanged = true ;
    end


% GET THE LIST OF AVAILABLE INPUTS
    function getAvailableSeedClasses()
        % Initialize
            availableSeedClasses = [] ;
        % Get navDIC_seeds folder path
            path = [hd.RootPath,'/navDIC_seeds'] ;
        % Get subfolders
            dirinfo = dir(path) ;
            dirinfo(~[dirinfo.isdir]) = [] ;
            dirinfo(ismember({dirinfo.name},{'.','..'})) = [] ;
        % Get Sub-classes
            for c = 1:length(dirinfo)
                seedClass = dirinfo(c).name ;
                % SubClasses dir
                    subClasses = dir([path,'/',seedClass]) ;
                    subClasses(~[subClasses.isdir]) = [] ;
                    subClasses(ismember({subClasses.name},{'.','..'})) = [] ;
                % Available class
                    superClass = regexprep(seedClass,'navDIC_','') ;
                    for sub = 1:length(subClasses)
                        subClass = regexprep(subClasses(sub).name,'navDIC_','') ;
                        availableSeedClasses(end+1).Name = [regexprep(superClass,'_',' | '),' | ',subClass] ;
                        availableSeedClasses(end).SuperClass = superClass ;
                        availableSeedClasses(end).Class = ['navDICSeed_',subClass] ;
                        availableSeedClasses(end).ClassName = subClass ;
                        availableSeedClasses(end).Folder = [path,'/navDIC_seeds/',seedClass,'/',subClasses(sub).name] ;
                    end
            end
        % Display Results
            disp('Available Seed Classes :') ;
            disp({availableSeedClasses.Name}') ;
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
            answer = questdlg(['ARE YOU SURE YOU WANT TO CLEAR ALL SEEDS ?']) ;
            if ~strcmp(upper(answer),'YES') ; return ; end
        % Clear all seeds
            resetSeeds() ;
        % Update Lists
            updateSeedList() ;
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
            fig = figure('tag','manageSeeds',...
                            'NumberTitle','off',...
                            'Name','MANAGE SEEDS',...
                            ...'windowstyle','modal',...
                            'toolbar','none',...
                            'menubar','none') ;
            fig.Units = 'normalized' ;
            fig.Position = [.5-figSize/2 .5-figSize/2 figSize figSize] ;
        % Second list
            listBoxSeeds = uicontrol(fig,'style','listbox') ;
            listBoxSeeds.Units = 'normalized' ;
            listBoxSeeds.Position = [margin 2*margin+btnHeight .5-2*margin-btnHeight/2 1-4*margin-btnHeight-ttlHeight] ;
        % First list Title
            ttl1 = uicontrol(fig,'style','text','string','DEFINED SEEDS','BackgroundColor','w','fontweight','bold','units','normalized') ;
            ttl1.Position = [(listBoxSeeds.Position(1)+sum(listBoxSeeds.Position([1,3])))/2-ttlWidth/2 1-margin-ttlHeight ttlWidth ttlHeight] ;
        % PopupMenu for choosing a subclass
            popupClasses = uicontrol(fig,'style','popupmenu') ;
            popupClasses.String = {availableSeedClasses.Name} ;
            popupClasses.Units = 'normalized' ;
            popupClasses.Position = [.5+1*margin+btnHeight/2 1-margin-ttlHeight .5-2*margin-btnHeight/2 ttlHeight] ;
        % Second list Title
            infoText = uicontrol(fig,'style','text','string','INFOS','BackgroundColor','w','fontweight','bold','units','normalized') ;
            infoText.Position = [.5+1*margin+btnHeight/2 2*margin+btnHeight .5-2*margin-btnHeight/2 1-4*margin-btnHeight-ttlHeight] ;
            infoText.HorizontalAlignment = 'left' ;
        % AddSeed Button
            btnAdd = uicontrol('style','pushbutton','units','normalized','string','+','fontweight','bold') ;
            btnAdd.Position = [.5-btnHeight/2 .5+btnHeight/2+margin/2 btnHeight btnHeight] ;
            btnAdd.Callback = @(src,evt)addSeed() ;
        % RemoveSeed Button
            btnRemove = uicontrol('style','pushbutton','units','normalized','string','-','fontweight','bold') ;
            btnRemove.Position = [.5-btnHeight/2 .5-btnHeight/2 btnHeight btnHeight] ;
            btnRemove.Callback = @(src,evt)removeSeed() ;
        % ModifySeed Button
            btnInfos = uicontrol('style','pushbutton','units','normalized','string','Modify','fontweight','bold') ;
            btnInfos.Position = [.5-btnHeight/2 .5-3*btnHeight/2-margin btnHeight btnHeight] ;
            btnInfos.Callback = @(src,evt)modifySeed() ;
        % DuplicateSeed Button
            btnInfos = uicontrol('style','pushbutton','units','normalized','string','Copy','fontweight','bold') ;
            btnInfos.Position = [.5-btnHeight/2 .5-5*btnHeight/2-2*margin btnHeight btnHeight] ;
            btnInfos.Callback = @(src,evt)duplicateSeed() ;
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
