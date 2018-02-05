
% POINTS
    file = fopen('Contour.txt') ;
    Vertices = str2num(fscanf(file,'%s')) ;
    Vertices = Vertices(:,1:2)/1000 ; % 2D --> 3D, et passage en mètres
    X = Vertices(:,1) ;
    Y = Vertices(:,2) ;
    nPts = length(X) ;
    fclose(file) ;  
    
clf
plot(X,Y,'k')
axis off
axis equal
    %%
pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
    1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];

decimP = 10 ;
pv = [X(1:decimP:end),Y(1:decimP:end)] ;
xmin = min(X) ;
xmax = max(X) ;
ymin = min(Y) ;
ymax = max(Y) ;
bbox = [xmin,ymin; xmax,ymax] ;
h0 = 0.02 ;


fstats=@(p,t) fprintf('%d nodes, %d elements, min quality %.2f\n', ...
                      size(p,1),size(t,1),min(simpqual(p,t)));

fprintf('Polygon\n');
[p,t]=distmesh2d(@dpoly,@huniform,h0,bbox,[],pv);
fstats(p,t);

