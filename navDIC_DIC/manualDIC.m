%% IMPLEMENT A MANUAL DIC

    % CLEAN WORKSPACE
        clc
        global hd
        clearvars -except hd

    % INITIALIZATION PARAMETERS
        camID = 1 ;
        seedNumber = 2 ;
        frames = '[1:end]' ; '[1,202,282]' ; '[20 42 43 72 73 78 79 84 85 126 127 168 169 186]' ; % Frames taken for DIC (allows decimation)
        dicDir = 1 ; % DIC running direction ('forward=1' or 'backward=-1')
        refFrame = 'first' ; % Reference image ('first' , 'last' or number)
        compConfig = 'Previous' ; % background help configuration ('Reference' or 'Previous') ;
        startWithNavDICPositions = 'all' ;
        addPreviousCorrection = false ; % When possible, add the previous correction (velocity or difference with navDIC positions) to the initialization
        exportTOnavDIC = false ;
        averagePreviousFrames = true ; % Ref frame is the average of the previous/next ones in forward/backward modes
        normToImageClassRange = true ; % Normalize images to their dataclass range
        timeMeanLength = 0 ; % Time averaging of images
        showInit = false ;
        codeProfile = false ; 
        figTag = 'Manual DIC' ;
        
    % LOAD AND PROCESS FRAMES
        globalDIC_01_LoadFrames ;
    % LOAD THE SEED
        Nodes = hd.Seeds(seedNumber).Points ; nNodes = size(Nodes,1) ;
        Elems = hd.Seeds(seedNumber).Triangles ;
        
    % PLOT THE REFERENCE CONFIGURATION
        btnHeight = 0.05 ;
        btnWidth = 0.05 ;
        margin = 0.003 ;
        clf(figDIC) ;
        ax = gobjects(0) ;
        ax(1) = mysubplot(1,2,1) ;
            im0 = imagesc(IMG(:,:,:,refFrame),'hittest','off') ;
            ax(1).XTick = [] ; ax(1).YTick = [] ;
            ax(1).YDir = 'reverse' ;
            box on
            axis equal
            axis tight 
            colormap gray
            title([compConfig,' Configuration'])
            mesh0 = patch('Faces',Elems,'Vertices',Nodes,'edgecolor','r','facecolor','none') ;
            pts0 = plot(ax(1),Nodes(:,1),Nodes(:,2),'.r','markersize',20) ;
        ax(2) = mysubplot(1,2,2) ;
            im = imagesc(IMG(:,:,:,refFrame),'hittest','off') ;
            ax(2).XTick = [] ; ax(2).YTick = [] ;
            ax(2).YDir = 'reverse' ;
            %ax(2).Clipping = 'off' ;
            box on
            axis equal
            axis tight 
            colormap gray
            ttl = title('Current Configuration') ;
            mesh = patch('Faces',Elems,'Vertices',Nodes,'edgecolor','b','facecolor','none') ;
            pts = impoint.empty ;
            for n = 1:size(Nodes,1)
                pts(end+1) = impoint(ax(2),Nodes(n,:)) ;
            end
            getPositions = @()reshape(cell2mat(arrayfun(@getPosition,pts,'UniformOutput',false)),2,[])' ;
            for n = 1:size(Nodes,1)
                addNewPositionCallback(pts(n),@(src,evt)set(mesh,'Vertices',getPositions())) ;
            end
        linkaxes(ax,'xy') ;
        stopBtn = uicontrol(figDIC,'style','togglebutton'...
                            ,'string','STOP'...
                            ,'units','normalized'...
                            ,'position',[1-btnWidth-margin margin btnWidth btnHeight]...
                            ,'callback',@(src,evt)disp('Stop!')) ;
        nextBtn = uicontrol(figDIC,'style','togglebutton'...
                        ,'string','NEXT'...
                        ,'units','normalized'...
                        ,'position',[1-2*btnWidth-2*margin margin btnWidth btnHeight]...
                        ,'callback',@(src,evt)disp('Next!')) ;
                    
    % MANUAL DIC LOOP
        % Motion Initialization
            useNavDICXn = false(nFrames,1) ;
            switch startWithNavDICPositions
                case 'all'
                    useNavDICXn(:) = true ;
                case 'none'
                    useNavDICXn(:) = false ;
                otherwise
                    useNavDICXn(startWithNavDICPositions) = true ;
            end
            Xn0 = ones([nNodes,2,nFrames])*NaN ;
            Xn0(:,:,useNavDICXn) = hd.Seeds(seedNumber).MovingPoints(:,:,frames(useNavDICXn)) ;
            Xn0(:,:,avgFrames) = repmat(Nodes,[1 1 length(avgFrames)]) ;
            Xn = Xn0 ; % Xn0 keeps a copy of the initialization
        % Loop
            for ii = dicFrames
                % First Guess for the positions
                    % Take the previously computed frame
                        if ~useNavDICXn(ii) || all(isnan(Xn(:,1,ii)))
                            Xn(:,:,ii) = Xn(:,:,ii-dicDir) ;
                        end
                    % Add the previous "correction" as convergence help
                        if addPreviousCorrection && abs(refFrame-ii)>=2
                            if useNavDICXn(ii-dicDir) % Add the correction of the previous frame with regard to the navDIC positions
                                correctionXn = Xn(:,:,ii-dicDir) - Xn0(:,:,ii-dicDir) ;
                            else % Add the correction of the previous frame with regard to the before-the-previous frame
                                correctionXn = (Xn(:,:,ii-dicDir)-Xn(:,:,ii-2*dicDir)) * (frames(ii)-frames(ii-dicDir))/(frames(ii-dicDir)-frames(ii-2*dicDir)) ;
                            end
                            Xn(:,:,ii) = Xn(:,:,ii) + correctionXn * 1.0 ;
                        end
                % Displacement guess
                    Un(:,:,ii) = Xn(:,:,ii) - Nodes ;
                % Update the config
                    im.CData = IMG(:,:,:,ii) ;
                    ttl.String = ['Adjust Current Configuration, Frame ',num2str(frames(ii))] ;
                    if strcmp(compConfig,'Previous')
                        im0.CData = IMG(:,:,:,ii-dicDir) ;
                        mesh0.Vertices = Xn(:,:,ii-dicDir) ;
                        pts0.XData = Xn(:,1,ii-dicDir) ;
                        pts0.YData = Xn(:,2,ii-dicDir) ;
                    end
                    %mesh.Vertices = Xn(:,:,ii) ;
                    for n = 1:size(Nodes,1)
                        setPosition(pts(n),Xn(n,:,ii)) ;
                    end
                % Let the user tune the positions
                    nextBtn.Value = false ;
                    while ~stopBtn.Value && ~nextBtn.Value
                        drawnow ;
                    end
                % Record the positions
                Xn(:,:,ii) = getPositions() ;
                if stopBtn.Value ; break ; end
            end
            ttl.String = 'Final Configuration' ;
            delete(nextBtn) ;
            drawnow ;
            
    % Post-process
        Un = Xn-Nodes ;
            
    % SEND THE RESULT TO navDIC
        if exportTOnavDIC

            if stopBtn.Value
                answer = questdlg({'PROCESS STOPPED BEFORE COMPLETION.', 'EXPORT THE RESULTS TO navDIC ?'},'/!\ WARNING','Yes','No','No') ;
                if isempty(answer) || strcmp(answer,'No') ; return ; end
            end

            if nNodes==1
                hd.Seeds(seedNumber).MovingPoints = permute(interp1(frames(:),permute(Xn,[3 2 1]),navDICFrames(:),'linear',NaN),[3 2 1]) ;
            else
            hd.Seeds(seedNumber).MovingPoints = interpn(...
                                                    repmat((1:nNodes)',[1 2 nFrames]),...
                                                    repmat(1:2,[nNodes 1 nFrames]),...
                                                    repmat(reshape(frames,[1 1 nFrames]),[nNodes 2 1]),...
                                                    Xn,...
                                                    repmat((1:nNodes)',[1 2 hd.nFrames]),...
                                                    repmat(1:2,[nNodes 1 hd.nFrames]),...
                                                    repmat(reshape(navDICFrames,[1 1 hd.nFrames]),[nNodes 2 1]),...
                                                'linear',NaN) ;
%             hd.Seeds(seedNumber).Displacements = interpn(...
%                                                     repmat((1:nNodes)',[1 2 nFrames]),...
%                                                     repmat(1:2,[nNodes 1 nFrames]),...
%                                                     repmat(reshape(frames,[1 1 nFrames]),[nNodes 2 1]),...
%                                                     Un,...
%                                                     repmat((1:nNodes)',[1 2 hd.nFrames]),...
%                                                     repmat(1:2,[nNodes 1 hd.nFrames]),...
%                                                     repmat(reshape(navDICFrames,[1 1 hd.nFrames]),[nNodes 2 1]),...
%                                                 'linear',NaN) ;
            end
        end
        

            
            
            
            
            
        
