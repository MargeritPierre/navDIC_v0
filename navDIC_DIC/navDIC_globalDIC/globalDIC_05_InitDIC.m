%% INITIALIZE DIC RESULTS

% INIT FIGURE
    infosHeight = 0.03 ;
    btnWidth = 0.05 ;
    margin = 0.003 ;
    fontSize = 12 ;
    figure(figGlobalDIC) ;
    figGlobalDIC = clf(figGlobalDIC,'reset') ;
    figGlobalDIC.Tag = 'globalDICfigure' ;
        ax = gobjects(0) ;
        ax(1) = mysubplot((nI<nJ)+1,(nI>=nJ)+1,1) ;
            ax(1).Position = ax(1).Position.*[1 1-infosHeight 1 1-infosHeight] ;
            im = imagesc(1:nJ,1:nI,Smooth(img0)) ;
            mesh = trisurf(Elems,Nodes(:,1),Nodes(:,2),Nodes(:,1)*0,'facecolor','none','edgecolor','r','linewidth',0.5,'edgealpha',0.5,'facealpha',0.5) ;
            markers = plot(NaN,NaN,'.b','markersize',15) ; % Deugging...
            colormap(ax(1),jet)
            set(ax(1),'Clipping','off') ;
        ax(2) = mysubplot((nI<nJ)+1,(nI>=nJ)+1,2) ;
            ax(2).Position = ax(2).Position.*[1 1-infosHeight 1 1-infosHeight] ;
            imRes = imagesc(1:nJ,1:nI,Smooth(img0)) ; 
            colormap(ax(2),gray)
            %ttl = title('','interpreter','none','units','normalized','color',[1 0 0]*1.0,'position',[0 1],'horizontalalignment','left','verticalalignment','top') ;
        axis(ax,'tight')
        axis(ax,'equal')
        set(ax,'xtick',[],'ytick',[])
        set(ax,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
        set(ax,'ydir','reverse') 
        clrbr = colorbar(ax(2)) ;
            clrbr.Location = 'west' ;
            clrbr.Color = 'w' ;
    infosText = uicontrol(figGlobalDIC,'style','text'...
                        ,'string','Infos: '...
                        ,'units','normalized'...
                        ,'backgroundcolor','w'...
                        ,'foregroundcolor','r'...
                        ,'fontsize',fontSize ...
                        ,'fontweight','bold' ...
                        ,'horizontalalignment','left'...
                        ,'position',[0 1-infosHeight 1 infosHeight]...
                        ) ;
    stopBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','STOP'...
                        ,'fontsize',fontSize/1.3 ...
                        ,'units','normalized'...
                        ,'position',[1-btnWidth-margin 1-infosHeight btnWidth infosHeight]...
                        ,'callback',@(src,evt)disp('Stop!')) ;
    if pauseAtPlot
        nextBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','NEXT'...
                        ,'units','normalized'...
                        ,'position',[1-2*btnWidth-2*margin 1-infosHeight btnWidth infosHeight]...
                        ,'callback',@(src,evt)disp('Next!')) ;
        continueBtn = uicontrol(figGlobalDIC,'style','togglebutton'...
                        ,'string','CONTINUE'...
                        ,'units','normalized'...
                        ,'position',[1-3*btnWidth-3*margin 1-infosHeight btnWidth infosHeight]...
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
    % Reference Image
        img1 = Smooth(img0) ;
        refImageChanged = true ;
    % Valid geometry masks
        VALID = true(nNodes,1) ;
        validElems = true(nElems,1) ;
        validEdges = true(nEdges,1) ;
        nakedEdges = sum(tri2edg(:,validElems),2)<2  ;
    