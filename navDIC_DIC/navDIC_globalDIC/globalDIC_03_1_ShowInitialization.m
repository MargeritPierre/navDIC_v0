  
%% SHOW THE INITIALIZATION THAT HAS BEEN ACHIEVED

        clf(figDIC) ;
        figure(figDIC) ;
        
        axes('nextplot','add','position',[0 0 1 1])
            im = imagesc(1:nJ,1:nI,img0) ; colormap(gray)
            ttl = title('','interpreter','none','units','normalized','position',[.005 0.995],'verticalalignment','top','horizontalalignment','left','color','r') ;
            mesh = patch('faces',Elems,'vertices',Nodes,'facecolor','none','edgecolor','r','linewidth',0.5,'edgealpha',0.5,'facealpha',0.5) ;
            axis tight
            axis equal
            set(gca,'xtick',[],'ytick',[])
            set(gca,'xlim',[0 nJ]+.5,'ylim',[0 nI]+.5)
            set(gca,'ydir','reverse')
            box on
            
        for ii = 1:nFrames
            im.CData = IMG(:,:,:,ii) ;
            ttl.String = num2str(ii) ;
            drawnow
            %pause(0.01)
        end
    