%% PROCESS THE SLICE IMAGES
clc
    
% RETRIEVE THE SLICE IMAGE
    sliceImg = findobj(groot,'tag','sliceImg') ;
    if isempty(sliceImg) ; warning('NO SLICE IMAGE FOUND') ; return ; end
    IMG0 = sliceImg.CData ;
    [nPts,nFrames,nColors] = size(IMG0) ;

% INITIALIZE A FIGURE
    % Is the figure already existing ?
        figTag = 'sliceProcessFigure' ;
            fig = findobj(groot,'tag',figTag) ;
                if isempty(fig)
                    fig = figure ;
                    fig.Position = [fig.Position(1:2)+fig.Position(3:4)/2 0 0] ...
                                    + 0.7*[-fig.Position(3:4)/2 fig.Position(3:4)] ;
                end
    % Reset the figure
        fig = clf(fig,'reset') ;
            fig.Tag = figTag ;
            fig.NumberTitle = 'off' ;
            fig.Name = 'Process Slice Figure' ;
        axImg = axes() ;
            axImg.YDir = 'reverse' ;
            axImg.LooseInset = 0.002*[1 1 1 1] ;
            axImg.XTick = [] ;
            axImg.YTick = [] ;
            axImg.XLim = [1 nFrames] ;
            axImg.YTick = [] ;
            axis off
            axis tight
        imgH = imagesc(IMG0) ;
            imgH.AlphaData = 1 ;
            colormap(axImg,'gray') ;
            
            
%% SELECT THE ROI

    imgH.CData = IMG0 ;
    imgH.AlphaData = 1 ;
    polyTag = 'maskPolygon' ;
    if exist('polygon') ; delete(polygon) ; end
    fcn = makeConstrainToRectFcn('impoly',[0 nFrames+1],[0 nPts+1]) ;
    polygon = impoly(axImg,'PositionConstraintFcn',fcn) ;
            
    
%% PROCESSING
    
    % Parameters
        % Pre-process
            applyROI = true ;
            setEndsToZero = false ;
            normalize = false ;
            extendEnds = false ;
        % Filter
            % Spatial Filter
                applySpatialFilter = false ;
                lambdaMin = 2 ; %
                lambdaMax = 2*10 ; %
                order = 5*10 ;
            % Time Filter (median filtering)
                applyTimeFilter = true ;
                    medFiltSize = 70 ;
        
    % Initialization
        IMG = double(IMG0) ;
        MASK = true(size(IMG)) ;
        
    % Pre-processing
        if applyROI
            if ~exist('polygon') || ~isvalid(polygon) ; warning('NO POLYGON DEFINING ROI') ; return ; end
            polyPos = getPosition(polygon) ;
            [JJ,II] = meshgrid(1:nFrames,1:nPts) ;
            MASK = inpolygon(JJ,II,polyPos(:,1),polyPos(:,2)) ;
            IMG(~MASK) = NaN ;
        end
        if setEndsToZero
            x = (0:nPts-1)'/(nPts-1) ;
            correct = (1-x)*IMG(1,:) + x*IMG(end,:) ;
            IMG = IMG - correct ;
        end
        if normalize
            IMG = bsxfun(@(x,y)x-y,IMG,median(IMG,1,'omitnan')) ; % Zero-mean
            %normCols = sqrt(sum(IMG.^2,1,'omitnan')./sum(MASK,1)) ;
            normCols = range(IMG,1) ;
            %IMG = bsxfun(@(x,y)x-y,IMG,normCols) ;
        end
        if extendEnds
            for pt = 1:nPts
                firstNotNaN = find(~isnan(IMG(pt,:)),1,'first') ;
                lastNotNaN = find(~isnan(IMG(pt,:)),1,'last') ;
                if ~isempty(firstNotNaN)
                    IMG(pt,1:firstNotNaN) = IMG(pt,firstNotNaN) ;
                    IMG(pt,lastNotNaN:end) = IMG(pt,lastNotNaN) ;
                end
            end
        end
        
    % Filtering
        if applySpatialFilter
            IMG(~MASK) = 0 ;
            filt = fir1(order,1./[lambdaMax lambdaMin],'bandpass') ;
            IMG = conv2(IMG,filt(:),'same') ;
            IMG(~MASK) = NaN ;
        end
        if applyTimeFilter
            IMG = medfilt1(IMG,medFiltSize,[],2,'omitnan','truncate') ;
        end
        
    % Display the result    
        imgH.CData = IMG ;
        imgH.AlphaData = MASK ;
        
        
%% PREPARE THE CORRELATION POINTS

    % PARAMETERS
        kern = blackman(9) ; % pre-filtering kernel
        
    % INITIALIZE A NEW FIGURE
        % Is the figure already existing ?
            figTag2 = 'sliceDICFigure' ;
                fig2 = findobj(groot,'tag',figTag2) ;
                    if isempty(fig2)
                        fig2 = figure ;
                        fig2.Position = [fig2.Position(1:2)+fig2.Position(3:4)/2 0 0] ...
                                        + 0.7*[-fig2.Position(3:4)/2 fig2.Position(3:4)] ;
                    end
        % Reset the figure
            figure(fig2) ;
                clf(fig2,'reset') ;
                fig2.Tag = figTag2 ;
                fig2.NumberTitle = 'off' ;
                fig2.Name = 'DIC on Slice' ;

    % Get the last full valid time (maximum of valid points)
        validPts = sum(MASK,1) ;
        maxValidPts = max(validPts) ;
        lastValidLine = find(validPts==maxValidPts,1,'last') ;
        
    % Plot the pixel intensities on this line
        lastLineValues = IMG(:,lastValidLine) ;
        filtVal = conv(lastLineValues,kern/sum(kern),'same') ;
        plot(lastLineValues) ;
        plot(filtVal) ;
        
    % Detect interest points (local minimas)
        d1 = (filtVal(3:end) - filtVal(1:end-2))/2 ;
        d2 = (filtVal(3:end) - 2*filtVal(2:end-1) + filtVal(1:end-2))/2 ;
        Nodes = find(abs(d1(2:end-1))<abs(d1(1:end-2)) & abs(d1(2:end-1))<abs(d1(3:end)) & d2(2:end-1)>0)+2 ;
        plot(Nodes,filtVal(Nodes),'.r','markersize',15)
        
    % Display informations
        nNodes = numel(Nodes) ;
        disp(newline)
        disp(['Number of Nodes Identified: ',num2str(nNodes)])
        disp(['Points Distance: ',num2str(median(diff(Nodes)))...
                ,'pix [',num2str(min(diff(Nodes)))...
                ,';',num2str(max(diff(Nodes)))...
                ,']' ...
                ])
            
%% PERFORM DIC

    % PARAMETERS
        clc
        plotEachIteration = false ;
        pauseExecution = false ;
        minNorm = 1e-5 ;
        maxIt = 500 ;
        minCorrCoeff = 0.1 ;
        
    % CONSTANT OBJECTS
        % Image
            dicIMG = IMG ;
            dicIMG(isnan(dicIMG)) = 0 ;
        % Shape Functions
            nElems = nNodes-1 ;
            MAPPING = zeros(nPts,nNodes) ;
            INSIDE = zeros(nPts,nElems) ;
            for elmt = 1:nElems
                pt1 = Nodes(elmt) ; pt2 = Nodes(elmt+1) ;
                x = (0:pt2-pt1)/(pt2-pt1) ;
                MAPPING(pt1:pt2,elmt) = 1-x ;
                MAPPING(pt1:pt2,elmt+1) = x ;
                INSIDE(pt1:pt2,elmt) = 1 ;
            end
            MAPPING = sparse(MAPPING) ;
            INSIDE = sparse(INSIDE) ;
            %figure(fig2) ; clf('reset') ; plot(MAPPING) ;
        % Reference line and derivatives
            F = dicIMG(:,lastValidLine) ;
            dF = conv(F,[1 0 -1]/2,'same') ;
            dF_da = sparse(dF(:).*MAPPING) ;
            Hess = dF_da'*dF_da ;
            
            
    % INITIALIZATION
        % Graphical objects
            figure(fig) ;
            fig.WindowButtonMotionFcn = '' ;
            if exist('dicPlot') ; delete(dicPlot) ; end
            dicPlot = plot(axImg,NaN,NaN,'r','linewidth',0.5) ;
            if exist('stopBtn') ; delete(stopBtn) ; end
            stopBtn = uicontrol(fig,'style','togglebutton'...
                                ,'string','STOP'...
                                ,'units','normalized'...
                                ,'position',[0.01 0.01 .08 .05]...
                                ,'callback',@(src,evt)disp('stop!')) ;
            if pauseExecution
                if exist('nextBtn') ; delete(nextBtn) ; end
                nextBtn = uicontrol(fig,'style','togglebutton'...
                                    ,'string','NEXT'...
                                    ,'units','normalized'...
                                    ,'position',[0.10 0.01 .08 .05]...
                                    ,'callback',@(src,evt)disp('next!')) ;
            end
            figure(fig2)
                clf ;
                plLin1 = plot(1:nPts,NaN*(1:nPts)) ;
                plLin2w = plot(1:nPts,NaN*(1:nPts)) ;
                plNodes = plot(Nodes,Nodes*NaN,'o') ;
                ttl = title('','interpreter','none','units','normalized','color',[1 0 0]*1.0,'position',[0 1],'horizontalalignment','left','verticalalignment','top') ;
        % Nodes Positions
            Xn = ones([nNodes,nFrames])*NaN ;
            Xn(:,lastValidLine) = Nodes ;
        % Displacements
            % Of Nodes
                Un = ones([nNodes,nFrames])*NaN ;
                Un(:,lastValidLine) = 0 ;
            % Of Pixels
                Up = zeros([nPts,1]) ;
        % Images
            lin1 = F ;
        % Mask
            VALID = true(nNodes,1) ;
            validElems = true(nElems,1) ;
                
    % RUN !
        for fr = lastValidLine-1:-1:1
            % Import and display image
                lin2 = dicIMG(:,fr) ;
            % Init
                it = 0 ;
                outFlag = false ;
                Un(:,fr) = Un(:,fr+1) ;
            % Add the previous "speed" as convergence help
                if fr<lastValidLine-1
                    Un(:,fr) = Un(:,fr) - (Un(:,fr+2)-Un(:,fr+1)) * 1.0 ;
                end
            % Newton-Raphson
                RMSE_0 = Inf ;
                lastPlotTime = tic ;
                while ~outFlag && ~stopBtn.Value
                    % Valid Domain
                        validElems = validElems & VALID(1:end-1) & VALID(2:end) ;
                        DOMAIN = any(INSIDE(:,validElems),2) ;
                    % IMAGE WARPING
                        % Compute the displacement at each pixel
                            Up = MAPPING(:,VALID)*Un(VALID,fr) ;
                        % Warp the image
                            lin2w = interp1(1:nPts,lin2,(1:nPts)+Up(:)','pchip',0) ;
                    % IMAGE MOMENTS
                        WEIGHT = INSIDE ; % MAPPING ; %
                        sumWEIGHT = sum(WEIGHT,1).' ;
                        % Mean over elements
                            meanImg1  = (WEIGHT'*lin1(:))./sumWEIGHT(:) ;
                            meanImg2w  = (WEIGHT'*lin2w(:))./sumWEIGHT(:) ;
                        % Zero-local-mean on pixels
                            lin1m = lin1(:)-WEIGHT*meanImg1(:) ;
                            lin2wm = lin2w(:)-WEIGHT*meanImg2w(:) ;
                        % Norm over element
                            normLin1 = sqrt(WEIGHT'*(lin1m(:).^2)) ;
                            normLin2w = sqrt(WEIGHT'*(lin2wm(:).^2)) ;
                        % Zero-local-mean-normalized images
                            lin1mz = lin1m(:)./(WEIGHT*normLin1) ;
                            lin2wmz = lin2wm(:)./(WEIGHT*normLin2w) ;
                    % IMAGE FUNCTIONAL
                        %diffImg = lin1(:)-lin2w(:) ; % Simple difference
                        %diffImg = lin1m(:)-lin2wm(:) ; % Zero-mean difference
                        diffImg = lin1mz(:)-lin2wmz(:) ; % Normalized Zero-mean difference
                    % NEWTON-RAPHSON PROCEDURE
                        % Compute the first RMSE derivative
                            dr_da = (diffImg(DOMAIN)'*dF_da(DOMAIN,VALID))' ;
                        % Updating DOFs
                            a = Hess(VALID,VALID)\dr_da ;
                        % Residues
                            %residueImg = norm(Hess(validDOF,validDOF)*a-dr_da) ;
                            %residueCONS = norm(CONS(validDOF,validDOF)*a+CONS(validDOF,validDOF)*[Un(VALID,1,ii);Un(VALID,2,ii)]) ;
                        % Displacement
                            Un(VALID,fr) = Un(VALID,fr) + a ;
                        % Positions
                            Xn(:,fr) = Nodes + Un(:,fr) ;
                    % VALIDATE GEOMETRY
                        % CULL OUT-OF-FRAME POINTS
                            VALID = VALID & Xn(:,fr)<=nPts & Xn(:,fr)>0 ;
                        % CULL OUT-OF-ROI points
                            VALID = VALID & MAPPING'*~MASK(:,fr)==0 ;
                        % Decorrelated elements
                            if minCorrCoeff>0
                                corrCoeff = abs(WEIGHT'*(lin1mz(:).*lin2wmz(:))) ;
                                switch size(WEIGHT,2) 
                                    case nElems % Correlation at the element level
                                        validElems = validElems & corrCoeff(:)>minCorrCoeff ; 
                                        VALID = VALID & any([[validElems ; 0] [0 ; validElems]],2) ;
                                    case nNodes % Correlation at the node level
                                        VALID = VALID & corrCoeff(:)>minCorrCoeff ;
                                end
                            end
                        % Set the non-valid values to NaN
                            Xn(~VALID,fr) = NaN ;
                            Un(~VALID,fr) = NaN ;
                    % CONVERGENCE CITERIONS
                        % Criterions
                            it = it+1 ;
                            %RMSE = norm(diffImg(DOMAIN))/norm(img1(DOMAIN)) ; 
                            normA = norm(a)/nNodes ; max(abs(a)) ;
                            nVALID = sum(VALID) ;
                        % Convergence criterion
                            if nVALID==0 ; outFlag = true ; end
                            if it>maxIt ; outFlag = true ; end
                            if normA<minNorm ; outFlag = true ; end
                            %if RMSE<1e-6 || abs((RMSE-RMSE_0)/RMSE) < 1e-4 ; outFlag = true ; end
                        % Keep the error
                            %RMSE_0 = RMSE ;
                    % DISPLAY
                        if plotEachIteration || (outFlag && toc(lastPlotTime)>0.05*0)
                            ttl.String = [num2str(fr),'(',num2str(it),')'] ;
                            dicPlot.XData = reshape(repmat([1:nFrames NaN]',[1 nNodes]),[],1) ;
                            dicPlot.YData = reshape([Xn' ; ones(1,nNodes)*NaN],[],1) ;
                            plLin1.YData = lin1mz ; plLin2w.YData = lin2wmz ;
                            %plLin1.YData = F ; plLin2w.YData = dF ;
                            %plLin1.YData = Up ; plNodes.YData = (Un(:,fr)-a).' ;
                            drawnow ;
                            lastPlotTime = tic ;
                        end
                    % Pause execution ?
                        if pauseExecution
                            while ~stopBtn.Value && ~nextBtn.Value
                                drawnow ;
                            end
                            nextBtn.Value = 0 ;
                        end
                end
                % Out Criterions
                    if stopBtn.Value; break ; end
                    if nVALID==0; break ; end
        end
        delete(stopBtn) ;
        if pauseExecution ; delete(nextBtn) ; end
        
    % New Reference Configuration
        lastValidFr = sum(isnan(Un(:,1:lastValidLine)),2)+1 ;
        lastValidPosition = Xn(sub2ind(size(Xn),(1:nNodes)',lastValidFr(:))) ;
        revUn = bsxfun(@minus,Xn,lastValidPosition) ;
        
        
        
%% RESULT INSPECTION
        
    uu = Un ;
    %uu = revUn ;
    
    figure(fig2) ;
        clf ;
        plDisp = plot(Nodes*NaN,Nodes,'o:') ;
        set(gca,'xlim',[min(uu(:)) max(uu(:))],'ylim',[1 nPts],'ydir','reverse')
        
    currFr = @()min(max(1,round(sum(sum(get(axImg,'CurrentPoint').*[1 0 0 ; 0 0 0])))),nFrames) ;
    showDisp = @()set(plDisp,'XData',uu(:,currFr())) ;
    
    fig.WindowButtonMotionFcn = @(src,evt)showDisp() ;
    
    
%% PLOT THE DISPLACEMENT
    
    fig.WindowButtonMotionFcn = '' ;
    
    figure(fig2)
        clf ;
        surf(revUn,'facecolor','interp') ;
        set(gca,'ydir','reverse')
        axis tight
        colorbar
        
        
    
        
    
    
    
    
    
    
    