clc
global hd

D = hd.Seeds(1).Displacements ;
D = reshape(D,[],hd.nFrames) ;
nanD = any(isnan(D),2) ;
D(nanD,:) = 0 ;
[U,S,V] = svd(D) ;

clf
plot(diag(S)) ;
set(gca,'yscale','log')

%%

P = hd.Seeds(1).Points ;
Tri = hd.Seeds(1).Triangles ;

r = 1 ;
amp = 0/100 ;

u = reshape(U(:,r),[],2) ;
u = u/max(abs(u(:))) ;
Lx = max(P(:,1))-min(P(:,1)) ;
Ly = max(P(:,2))-min(P(:,2)) ;
L = norm([Lx Ly]) ;

clf ;
ax(1) = mysubplot(1,2,1) ;
    trisurf(Tri,P(:,1)+u(:,1)*L*amp,P(:,2)+u(:,2)*L*amp,u(:,1),u(:,1))
    axis tight
    axis off
    myaxisequal('xy')
ax(2) = mysubplot(1,2,2) ;
    trisurf(Tri,P(:,1)+u(:,1)*L*amp,P(:,2)+u(:,2)*L*amp,u(:,2),u(:,2))
    axis tight
    axis off
    myaxisequal('xy')
    
linkprop(ax,'view') ;

