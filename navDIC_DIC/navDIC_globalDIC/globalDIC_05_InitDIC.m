%% INITIALIZE DIC RESULTS

% INIT FIGURE
    figure(figGlobalDIC) ;
    figGlobalDIC = clf(figGlobalDIC,'reset') ;
    figGlobalDIC.Tag = 'globalDICfigure' ;
        ax = [] ;
        ax(1) = mysubplot((nI<nJ)+1,(nI>=nJ)+1,1) ;
            im = imagesc(1:nJ,1:nI,Func(img0)) ;
            mesh = trisurf(Elems,Nodes(:,1),Nodes(:,2),Nodes(:,1)*0,'facecolor','none','edgecolor','r','linewidth',0.5,'edgealpha',0.5,'facealpha',0.5) ;
            markers = plot(NaN,NaN,'.b','markersize',15) ; % Deugging...
            colormap(ax(1),jet)
            set(ax(1),'Clipping','off') ;
        ax(2) = mysubplot((nI<nJ)+1,(nI>=nJ)+1,2) ;
            imRes = imagesc(1:nJ,1:nI,Func(img0)) ; 
            colormap(ax(2),gray)
            ttl = title('','interpreter','none','units','normalized','color',[1 0 0]*1.0,'position',[0 1],'horizontalalignment','left','verticalalignment','top') ;
        axis(ax,'tight')
        axis(ax,'equal')
        set(ax,'xtick',[],'ytick',[])
        set(ax,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
        set(ax,'ydir','reverse')
    infosText = uicontrol(figGlobalDIC,'style','text'...
                        ,'string','Infos: '...
                        ,'units','normalized'...
                        ,'backgroundcolor','w'...
                        ,'foregroundcolor','r'...
                        ,'fontsize',12 ...
                        ,'fontweight','bold' ...
                        ,'horizontalalignment','left'...
                        ,'position',[0 .95 1 .05]...
                        ) ;
    stopBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','STOP'...
                        ,'units','normalized'...
                        ,'position',[0.01 0.01 .08 .05]...
                        ,'callback',@(src,evt)disp('Stop!')) ;
    if pauseAtPlot
        nextBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','NEXT'...
                        ,'units','normalized'...
                        ,'position',[0.1 0.01 .08 .05]...
                        ,'callback',@(src,evt)disp('Next!')) ;
        continueBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','CONTINUE'...
                        ,'units','normalized'...
                        ,'position',[0.19 0.01 .08 .05]...
                        ,'callback',@(src,evt)disp('Continue!')) ;
    end
    
% Other figure if needed to debug
    if 0
        if isempty(findobj(0,'tag','figDebug'))
            figDebug = figure('tag','figDebug') ;
        end
        figDebug = figure(findobj(0,'tag','figDebug')) ;
    end
        

% INITIALIZE
    % Nodes position
        if startWithNavDICPositions
            Xn = hd.Seeds(seedNumber).MovingPoints(:,:,frames) ;
            Xn(:,:,refFrame) = Nodes ;
        else
            Xn = ones([nNodes,2,nFrames])*NaN ;
            Xn(:,:,avgFrames) = repmat(Nodes,[1 1 length(avgFrames)]) ;
        end
    % Displacements
        % Of Nodes
            Un = Xn - Nodes ;
        % Of Pixels
            Up = zeros([nI nJ 2]) ;
    % Reference ImageGeome
        img1 = F ;
        % Image Vector
            img1v = img1(:) ;
        % Moments
            % Integration weights
                sumWEIGHT = sum(localWEIGHT,1).' ;
            % Mean over elements
                meanImg1 = (localWEIGHT'*img1v)./sumWEIGHT(:) ;
            % Zero-local-mean on pixels
                img1m = img1v-localWEIGHT*meanImg1(:) ;
            % Norm over element
                normImg1 = sqrt(localWEIGHT'*(img1m(:).^2)) ;
            % Zero-local-mean-normalized images
                img1mz = img1m(:)./(localWEIGHT*normImg1) ;
    % Mask
        VALID = true(nNodes,1) ;
        validElems = true(nElems,1) ;
        validEdges = true(nEdges,1) ;
        nakedEdges = sum(tri2edg(:,validElems),2)<2  ;
    