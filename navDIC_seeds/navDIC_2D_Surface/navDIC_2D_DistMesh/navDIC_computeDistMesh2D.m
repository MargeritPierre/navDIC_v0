function mesh = navDIC_computeDistMesh2D(fd,fh,h0,bbox,pfix,Axes)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
%   Adapted for the use in navDIC by Pierre Margerit (2018)

% Parameters for convergence of the mesh
    maxCount = 200 ; % maxNumber of iterations
    dptol=.001; 
    ttol=.01; 
    Fscale=1.2; 
    deltat=.05; 
    geps=.001*h0; 
    deps=sqrt(eps)*h0;
    densityctrlfreq=30;

% 1. Create initial distribution in bounding box (equilateral triangles)
    [x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
    x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
    p=[x(:),y(:)];                                       % List of node coordinates

% 2. Remove points outside the region, apply the rejection method
    % 2.1 Remove points outside BBOX
        p = p(p(:,1)>=bbox(1,1)&p(:,2)>=bbox(1,2)&p(:,1)<=bbox(2,1)&p(:,2)<=bbox(2,2),:) ;
        p=p(fd(p)<geps,:);                                   % Keep only d<0 points
        r0=1./fh(p).^2;                                      % Probability to keep point
        p=p(rand(size(p,1),1)<r0./max(r0),:);                % Rejection method
        if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end     % Remove duplicated nodes
        pfix=unique(pfix,'rows'); nfix=size(pfix,1);
        p=[pfix; p];                                         % Prepend fix points
        N=size(p,1);                                         % Number of points N

% Initialize data
    mesh = [] ;
    lastPlotTime = tic ;
    count=0;
    pold=inf;                                           % For first iteration

% Initialize mesh
    triMesh = findobj(Axes,'tag','DistMeshPreview') ;
    delete(triMesh) ;
    p0 = p ;
    t0 = computeMesh(p0) ;
    triMesh = patch('vertices',p0...
                    ,'faces',t0...
                    ,'edgecolor','b'...
                    ,'facecolor',flip([.8,.9,1])...
                    ,'tag','DistMeshPreview'...
                    ,'facealpha',0.7 ...
                    ,'hittest','off' ...
                    ) ;
    uistack(triMesh,'bottom') ;
    uistack(findobj(Axes,'type','image'),'bottom') 
    
while 1
  count=count+1;
  % 3. Retriangulation by the Delaunay algorithm
      if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
        pold=p;                                          % Save current positions
        t = computeMesh(p) ;
        % 4. Describe each bar by a unique pair of nodes
            bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
            bars=unique(sort(bars,2),'rows');                % Bars as node pairs
        % 5. Graphical output of the current mesh
            updateMeshPlot(p,t) ;
      end

  % 6. Move mesh points based on bar lengths L and forces F
      barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
      L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
      hbars=fh((p(bars(:,1),:)+p(bars(:,2),:))/2) ;
      L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
  
  % Density control - remove points that are too close
      if mod(count,densityctrlfreq)==0 && any(L0>2*L)
          p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix),:)=[];
          N=size(p,1); pold=inf;
          continue;
      end
      
  % Perturbations computation
      F=max(L0-L,0);                                     % Bar forces (scalars)
      Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)
      Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
      Ftot(1:size(pfix,1),:)=0;                          % Force = 0 at fixed points
      p=p+deltat*Ftot;                                   % Update node positions

  % 7. Bring outside points back to the boundary
      d=fd(p) ; ix=d>0;                 % Find points outside (d>0)
      dgradx=(fd([p(ix,1)+deps,p(ix,2)])-d(ix))/deps; % Numerical
      dgrady=(fd([p(ix,1),p(ix,2)+deps])-d(ix))/deps; %    gradient
      dgrad2=dgradx.^2+dgrady.^2;
      p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];    % Project

  % 8. Termination criterion: All interior nodes move less than dptol (scaled)
      if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end
      if ~isvalid(triMesh) ; break ; end
      if count>maxCount ; break ; end 
end

% Clean up and plot final mesh
[p,t]=fixmesh(p,t);

% Final Iteration
    updateMeshPlot(p,t) ;
    triMesh.EdgeColor = 'b' ;
    triMesh.FaceColor = [.8,.9,1] ;
    drawnow ;
    

% Format the output mesh
    mesh.Points = p ;
    mesh.Triangles = t ;
    mesh.Patches = triMesh ;
    
    
    
%simpplot(p,t)
    

         
% ===================================================================================================================    
% CHILD FUNCTIONS
% ===================================================================================================================

   
    function t = computeMesh(p)
        t=delaunay(p);                                  % List of triangles
        pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
        t=t(fd(pmid)<-geps,:);                          % Keep interior triangles
    end

    function updateMeshPlot(p,t)
        if isvalid(triMesh) 
            triMesh.Vertices = p ;
            triMesh.Faces = t ;
            if toc(lastPlotTime)>0.02
                drawnow
                lastPlotTime = tic ;
            end
        end
    end

end
