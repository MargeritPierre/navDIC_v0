%% COMPUTES THE DISPLACEMENT FIELDS CORRESPONDING TO PROBLEM SENSITIVITIES TO MATERIAL PARAMETERS

    % CLEAN WORKSPACE
        clc
        global hd
        clearvars -except hd
        
    % PARAMETERS
        seedNumber = 1 ;
        
    % GET THE GEOMETRIES AND PROCESS
        % Retrieve
            Nodes = hd.Seeds(seedNumber).Points ;
            Elems = hd.Seeds(seedNumber).Triangles ;
        % Sort Element Nodes in Trigo Order
            ElemNodes = reshape(cat(3,Nodes(Elems,1),Nodes(Elems,2)),[],3,2) ;
            Centers = mean(ElemNodes,2) ;
            RelPos = ElemNodes-Centers ;
            theta = angle(RelPos(:,:,1)+1i*RelPos(:,:,2)) ;
            [theta,indSort] = sort(theta,2,'ascend') ;
            Elems = Elems(sub2ind(size(Elems),repmat((1:size(Elems,1))',[1 3]),indSort)) ;
        % Edges creation and finding
            [Edges,~,elem2Edge] = unique(sort(reshape(cat(3,Elems,circshift(Elems,1,2)),[],2),2),'rows') ;
            elem2Edge = reshape(elem2Edge,[],3) ;
            elem2Edge = sparse(elem2Edge,repmat([1:size(Elems,1)]',[1 3]),1) ;
            NakedEdges = find(sum(elem2Edge,2)==1) ;
            EdgeNodes = permute(reshape(Nodes(Edges,:),[],2,2),[1 2 3]) ;
        % Outgoing normals on Naked Edges
            NakedNodes = EdgeNodes(NakedEdges,:,:) ;
            Pe = mean(NakedNodes,2) ; % Naked Edges middle point
            Normals = cat(3,-diff(NakedNodes(:,:,2),1,2),diff(NakedNodes(:,:,1),1,2)) ;
            Normals = Normals./sqrt(sum(Normals.^2,3)) ;
        % Re-orient trought outside direction
            Pt = mean(reshape(Nodes(elem2Edge(NakedEdges,:)*Elems,:),[],3,2),2) ;
            Normals = Normals.*sign(sum((Pe-Pt).*Normals,3)) ;
        % Build the operator (one eq. by free edge)
            Nx = sparse(repmat((1:length(NakedEdges))',[1 2]),Edges(NakedEdges,:),Normals(:,:,1)*[.5 .5]) ;
            Ny = sparse(repmat((1:length(NakedEdges))',[1 2]),Edges(NakedEdges,:),Normals(:,:,2)*[.5 .5]) ;
            
        
    % GRADIENT OPERATOR (Dx*f = df_dx; Dy*f = df_dy)
        Xt = reshape(Nodes(Elems',1),3,[])' ; Yt = reshape(Nodes(Elems',2),3,[])' ;
        areas = 1/2*(Xt(:,2).*Yt(:,3) - Xt(:,3).*Yt(:,2) - Xt(:,1).*Yt(:,3) + Xt(:,3).*Yt(:,1) + Xt(:,1).*Yt(:,2) - Xt(:,2).*Yt(:,1)) ;
        aa = cross(Xt,Yt,2) ;
        bb = circshift(Yt,2,2) - circshift(Yt,1,2) ;
        cc = circshift(Xt,1,2) - circshift(Xt,2,2) ;
        Dx = sparse(repmat((1:size(Elems,1))',[1 3]),Elems,1./2./areas(:).*bb,size(Elems,1),size(Nodes,1)) ;
        Dy = sparse(repmat((1:size(Elems,1))',[1 3]),Elems,1./2./areas(:).*cc,size(Elems,1),size(Nodes,1)) ;
        
    % INIT FIGURE
        figTag = 'StressModes' ;
        fig = findobj(0,'tag',figTag) ;
        fig = clf(fig) ;
        mesh = patch('Faces',Elems,'Vertices',Nodes,'FaceVertexCData',Nodes(:,1)*NaN,'facecolor','interp','edgecolor','k') ;
        naked = plot(EdgeNodes(NakedEdges,:,1)',EdgeNodes(NakedEdges,:,2)','b') ;
        arrowsFree = quiver(Pe(:,1),Pe(:,2),Normals(:,:,1),Normals(:,:,2),'b') ;
        arrowsLoad = quiver(Pe(:,1)*NaN,Pe(:,2)*NaN,Normals(:,:,1)*NaN,Normals(:,:,2)*NaN,'r') ;
        set(gca,'xtick',[],'ytick',[])
        axis equal
        axis tight
        axis off
        
    % LET THE USER CHOOSE THE LOADED POINTS
        title('SELECT THE LOADED EDGES (Quit with Right-Click)')
        Xe = mean(EdgeNodes(NakedEdges,:,1),2) ; Ye = mean(EdgeNodes(NakedEdges,:,2),2) ; 
        FreeEdges = NakedEdges ;
        while 1
            [x,y,button] = ginput(1) ;
            if button~=1 ; break ; end
            dist = (x-Xe).^2 + (y-Ye).^2 ;
            indLoad = find(dist==min(dist)) ;
            FreeEdges = setdiff(FreeEdges,NakedEdges(indLoad)) ;
            naked(indLoad).Color = 'r' ;
            arrowsFree.XData(indLoad) = NaN ;
            arrowsLoad.XData(indLoad) = Pe(indLoad,1) ;
            arrowsLoad.YData(indLoad) = Pe(indLoad,2) ;
            arrowsLoad.UData(indLoad) = Normals(indLoad,:,1) ;
            arrowsLoad.VData(indLoad) = Normals(indLoad,:,2) ;
        end
        
        
%% BUILD STIFFNESS MATRICES

    % Global Scale
        pixByMM = (928.8-568.1)/7 ;

    % Strain Operator B (linearized Green-Lagrange)
        O = Dx*0 ;
        B = [Dx O ; O Dy ; Dy Dx] ;
        B = B([1:3:end-2,2:3:end-1,3:3:end],[1:2:end-1,2:2:end]) ;
    
    % Elementary stiffness matrices
        K11 = B'*kron(spdiags(areas,1,length(areas),length(areas)),[1 0 0 ; 0 0 0 ; 0 0 0])*B ;
        K12 = B'*kron(spdiags(areas,1,length(areas),length(areas)),[0 1 0 ; 1 0 0 ; 0 0 0])*B ;
        K16 = B'*kron(spdiags(areas,1,length(areas),length(areas)),[0 0 1 ; 0 0 0 ; 1 0 0])*B ;
        K22 = B'*kron(spdiags(areas,1,length(areas),length(areas)),[0 0 0 ; 0 1 0 ; 0 0 0])*B ;
        K26 = B'*kron(spdiags(areas,1,length(areas),length(areas)),[0 0 0 ; 0 0 1 ; 0 1 0])*B ;
        K66 = B'*kron(spdiags(areas,1,length(areas),length(areas)),[0 0 0 ; 0 0 0 ; 0 0 1])*B ;
        Ku = speye(size(Nodes,1)) ;
        
    
    
    
        
    %% PLOT THE RESULT
    
        mode = 1  ;
        
        clf(fig) ;
        ax = gobjects(0) ;
        meshV = gobjects(0) ;
        for c = 1:size(V,2) 
            ax(c) = mysubplot(size(V,2),1,c) ;
                meshV(c) = patch('Vertices',Nodes,'Faces',Elems,'facecolor','interp','facealpha',0.5,'edgecolor','k','edgealpha',0.5) ;
                meshV(c).FaceVertexCData = V(:,c,mode) ;
                meshV(c).Vertices(:,3) = V(:,c,mode) ;
                colorbar
                myaxisequal('xy')
        end
        set(ax,'xtick',[],'ytick',[])
        %axis(ax,'equal')
        axis(ax,'tight')
        axis(ax,'off')
        
    %% LOWER-NORM SOLUTION
        A = = [Adiv ; Afree ; Ator ; Afor ;]
        
        
        
        
        