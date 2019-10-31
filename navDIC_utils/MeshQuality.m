

global hd

qmin = 1/10 ;
seedNumber = 3
mesh = hd.Seeds(seedNumber) ;
tri = mesh.Triangles ;
Xm = mesh.Points(:,1) ;
Ym = mesh.Points(:,2) ;

    % Quality
        bars = [tri(:,[1,2]);tri(:,[2,3]);tri(:,[3,1])] ;
        dX = diff(Xm(bars),1,2) ;
        dY = diff(Ym(bars),1,2) ;
        barLengths = sqrt(dX.^2+dY.^2) ;
        nTri = size(tri,1) ;
        a = barLengths(1:nTri) ;
        b = barLengths(nTri+(1:nTri)) ;
        c = barLengths(2*nTri+(1:nTri)) ;
        q = (b+c-a).*(c+a-b).*(a+b-c)./(a.*b.*c) ;
        
clf
axis equal
patch('Vertices',mesh.Points,'Faces',mesh.Triangles,'facecolor','none','edgecolor','k')
patch('Vertices',mesh.Points,'Faces',mesh.Triangles(q<qmin,:),'facecolor','r','edgecolor','r')

disp([num2str(sum(q<qmin)) ' elements with a quality less than ' num2str(qmin)]) ;

%% DELETE TRIANGLES

hd.Seeds(seedNumber).Triangles(q<qmin,:) = [] ;