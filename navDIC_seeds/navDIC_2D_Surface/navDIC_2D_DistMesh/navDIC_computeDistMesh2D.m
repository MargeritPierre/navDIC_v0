function mesh = navDIC_computeDistMesh2D(fd,fh,h0,bbox,pfix,Axes)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
%   Adapted for the use in navDIC by Pierre Margerit (2018)

% Parameters for convergence of the mesh
    maxCount = 1000 ; % max number of iterations
    dptol = .001; % point displacement tolerance (convergence criterion, relative to mean density)
    ttol = .2; % re-triangulation tolerance
    Fscale = 1.2; % spring L0 relative to bar length
    deltat = .3; % explicite time increment
    geps = .01*h0; % maximum distance function allowed (remove pts)
    deps = sqrt(eps)*h0; % discrete gradient estim. step
    densityctrlfreq = Inf ; % control if points too close every (.) iterations
    tooCloseThrs = 0.75 ;
    plotFreq = 10 ; % update plot frequency

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
    pold=-inf; % For first iteration
    nReTri = 0 ;
    infos = {}; % to print informations at the end
    
% Delete previous graphic handles
    gHandles = findobj(Axes,'tag','DistMeshPreview') ;
    delete(gHandles) ;

% Initialize mesh
    p0 = p ;
    t0 = computeMesh(p0) ;
    triMesh = patch('vertices',p0...
                    ,'faces',t0...
                    ,'edgecolor',[1 1 1]*0 ...
                    ,'facecolor','interp'...
                    ,'tag','DistMeshPreview'...
                    ,'facealpha',0.5 ...
                    ,'hittest','off' ...
                    ) ;
    uistack(triMesh,'bottom') ;
    uistack(findobj(Axes,'type','image'),'bottom') 
    
% Superimposing Interaction Axes
    infosText = uicontrol(Axes.Parent,'style','Text'...
                        ,'tag','DistMeshPreview' ...
                        ,'string',''...
                        ,'backgroundcolor','w' ...
                        ,'foregroundcolor','b' ...
                        ,'fontsize',10 ...
                        ,'fontweight','bold' ...
                        ,'units','normalized'...
                        ,'position',[0.005 1-0.045 0.99 .04]...
                        ,'horizontalalignment','left' ...
                        ) ;
    stopBtn = uicontrol(Axes.Parent,'style','togglebutton'...
                        ,'tag','DistMeshPreview' ...
                        ,'string','STOP'...
                        ,'units','normalized'...
                        ,'position',[1-0.055 1-0.055 .05 .05]...
                        ,'callback',@(src,evt)disp('Stop!')) ;
    
while 1
    % Initialize
        dP = 0 ;
        p0 = p ;
        dp = p-p0 ;
        count=count+1;
        
  % Retriangulation by the Delaunay algorithm
      %if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
      if any(isinf(pold(:))) || any(abs(sqrt(sum((p-pold).^2,2))./fh(p))>ttol)          % Any large movement?
        nReTri = nReTri+1 ;
        pold=p;                                          % Save current positions
        t = computeMesh(p) ;
        % 4. Describe each bar by a unique pair of nodes
            bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
            bars=unique(sort(bars,2),'rows');                % Bars as node pairs
        % If only fixed points remains, break
            if size(p,1)<=size(pfix,1) ; infos{end+1} = {'out criterion: only fixed points'} ; break ; end
      end

  % Move mesh points based on bar lengths L and forces F
      barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
      L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
      Lt=fh((p(bars(:,1),:)+p(bars(:,2),:))/2) ;         % Lt = desired bar length
      L0=Lt*Fscale*sqrt(sum(L.^2)/sum(Lt.^2));           % L0 = target lengths for update
  
  % Density control - remove points that are too close
      tooClose = L./Lt<tooCloseThrs ;
      if mod(count,densityctrlfreq)==0 && any(tooClose)
          p(setdiff(reshape(bars(tooClose,:),[],1),1:nfix),:)=[];
          N=size(p,1); pold=inf;
          continue;
      end
      
  % Perturbations computation
      F=max(L0-L,0);                                     % Bar forces (scalars)
      Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)
      Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
      Ftot(1:size(pfix,1),:)=0;                          % Force = 0 at fixed points
      p=p+deltat*Ftot ;                                  % Update node positions

  % 7. Bring outside points back to the boundary
      d=fd(p) ; ix=d>0;                 % Find points outside (d>0)
      dgradx=(fd([p(ix,1)+deps,p(ix,2)])-d(ix))/deps; % Numerical
      dgrady=(fd([p(ix,1),p(ix,2)+deps])-d(ix))/deps; %    gradient
      dgrad2=dgradx.^2+dgrady.^2;
      p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];    % Project

  % 8. Termination criterion: All interior nodes move less than dptol
  % (scaled)navDIC
      dp = p-p0 ;
      %dP = max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0) ;
      dP = max(sqrt(sum(dp.^2,2))/h0) ;
      if dP<dptol ; infos{end+1} = {'out criterion: |dP|<tol'} ; break; end
      if ~isvalid(triMesh) ; infos{end+1} = {'out criterion: TriMesh not valid'} ; break ; end
      if count>=maxCount ; infos{end+1} = {'out criterion: iteration count'} ; break ; end 
      if stopBtn.Value ; infos{end+1} = {'out criterion: user stop'} ; break ; end 

  % 5. Graphical output of the current mesh
      if toc(lastPlotTime)>1/plotFreq
          updateMeshPlot(p,t) ;
          lastPlotTime = tic ;
      end
end
infos{end+1} = {[num2str(count),' iterations']} ;
infos{end+1} = {['last dP: ',num2str(dP)]} ;
infos{end+1} = {[num2str(nReTri),' re-triangulations']} ;

% Clean up and plot final mesh
[p,t]=fixmesh(p,t);

% Final Iteration
    updateMeshPlot(p,t) ;
    triMesh.EdgeColor = 'b' ;
    triMesh.FaceColor = [.8,.9,1] ;
    delete(infosText) ;
    delete(stopBtn) ;
    drawnow ;
    

% Format the output mesh
    mesh.Points = p ;
    mesh.Triangles = t ;
    mesh.Patches = triMesh ;
    infos{end+1} = {[num2str(size(p,1)),' nodes']} ;
    infos{end+1} = {[num2str(size(t,1)),' triangles']} ;
    infos{end+1} = '<a href="matlab: opentoline(which(''navDIC_computeDistMesh2D.m''),6)">method parameters</a>' ;
    
    
% DISPLAY INFOS
disp(newline)
disp('--- DistMesh ---')
for in = [infos{:}]
    disp(['  ',in{1}])
end
disp('----------------')
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
            triMesh.FaceVertexCData = sqrt(sum(dp.^2,2)) ;
            triMesh.Vertices = p ;
            triMesh.Faces = t ;
            infosText.String = ['Infos: ' ...
                                 , ' | it: ' , num2str(count),'/',num2str(maxCount) ...
                                 , ' | Nodes: ' , num2str(length(p)) ...
                                 , ' | Triangles: ' , num2str(size(t,1)) ...
                                 , ' | Re-Tri: ' , num2str(nReTri) ...
                                 , ' | dP: ' , num2str(dP,3),'/',num2str(dptol,3) ...
                                ] ; 
            drawnow
        end
    end

end
