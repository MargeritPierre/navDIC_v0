%% TRY BOX INTERACTIVITY
clc

clf ;
ax = axes() ;
    ax.XTick = [] ;
    ax.YTick = [] ;
    box(ax,'on') ;
    axis(ax,'equal') ;
    ax.XLim = [0 1] ;
    ax.YLim = [0 1] ;
    
Lx = .1 ; Ly = .05 ;
X0 = .5 ; Y0 = .5 ;
xx = [-1 1 1 -1]*Lx + X0 ; yy = [-1 -1 1 1]*Ly + Y0 ;
%pa = patch(xx,yy,xx*0,'r') ;

pos = [xx(1) yy(1) range(xx) range(yy)] ;
Boxes = imrect(ax,pos) ;

for b = 1:length(Boxes) ; setResizable(Boxes(b),false) ; end


