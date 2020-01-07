global hd

Seed = hd.Seeds(end) ;

Pts = Seed.Points ; Seed.MovingPoints(:,:,1) ;
IMG = repmat(Seed.refImgs{1},[1 1 3]) ;

lbls = arrayfun(@(i)num2str(i),1:size(Pts,1),'UniformOutput',false) ;

clf reset ;
    im = image(IMG) ;
        axis ij
        axis equal
        axis tight
    pl = plot(Pts(:,1),Pts(:,2),'.r','markersize',20) ;
    txt = text(Pts(:,1),Pts(:,2),lbls) ;
        set(txt,'Interpreter','tex') ;
        set(txt,'Color','r') ;
        set(txt,'HorizontalAlignment','right','VerticalAlignment','bottom') ;
        set(txt,'FontSize',17) ;
        
        
        
        
        
        