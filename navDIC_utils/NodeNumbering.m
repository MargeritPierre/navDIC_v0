global hd

Seed = hd.Seeds(end) ;
cam = 1 ;
frame = hd.CurrentFrame ; 2000 ;

Pts = Seed.Points ; Seed.MovingPoints(:,:,frame) ;
Elems = Seed.Elems ;
IMG = hd.Images{cam}{frame} ;

% 
clf reset ;
    im = imagesc(IMG) ;
        axis ij
        axis equal
        axis tight
        colormap gray
    % Elements Labelling
        nod2elem = Seed.elem2nod('ones')' ; nod2elem = nod2elem./sum(nod2elem,2) ;
        C = nod2elem*Pts ;
        lbls = arrayfun(@(i)num2str(i),1:size(Elems,1),'UniformOutput',false) ;
        pl = patch('Vertices',Pts,'Faces',Elems,'EdgeColor','b','FaceColor','w','FaceAlpha',0.5) ;
        txt = text(C(:,1),C(:,2),lbls) ;
            set(txt,'Interpreter','tex') ;
            set(txt,'Color','b') ;
            set(txt,'HorizontalAlignment','center','VerticalAlignment','middle') ;
            set(txt,'FontSize',17) ;
    % Node Labelling
        lbls = arrayfun(@(i)num2str(i),1:size(Pts,1),'UniformOutput',false) ;
        pl = plot(Pts(:,1),Pts(:,2),'.r','markersize',20) ;
        txt = text(Pts(:,1),Pts(:,2),lbls) ;
            set(txt,'Interpreter','tex') ;
            set(txt,'Color','r') ;
            set(txt,'HorizontalAlignment','right','VerticalAlignment','bottom') ;
            set(txt,'FontSize',17) ;
        
        
        
        
        
        