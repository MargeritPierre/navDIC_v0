%% PROCESSING ELEMENT GEOMETRY
    Xt = Nodes(:,1) ; Yt = Nodes(:,2) ; % Separate coordinates for indexing
    Xc = mean(Xt(Elems),2) ; Yc = mean(Yt(Elems),2) ; % Barycentric coordinates
    Areas = zeros(nElems,1) ; % Element areas
    %gradT = zeros(nElems,2*nNodes,2,2) ; % Gradient matrix (element-wise)
    tri2nod = zeros(nNodes,nElems) ; % triangles to nodes transform
    % Gradient values
        iiiD = zeros(nElems,3) ;
        jjjD = zeros(nElems,3) ;
        vvvD1 = zeros(nElems,3) ;
        vvvD2 = zeros(nElems,3) ;
    wtbr = waitbar(0,'Processing elements') ;
    ticTime = tic ;
    for tri = 1:nElems
        % Points associated to the element
            Xtri = Nodes(Elems(tri,:),1) ;
            Ytri = Nodes(Elems(tri,:),2) ;
        % Sort nodes trigonometric orientation
            pol = angle(Xtri-Xc(tri)) + 1i*(Ytri-Yc(tri)) ;
            [~,ind] = sort(pol,'ascend') ;
            Elems(tri,:) = Elems(tri,ind) ;
            Xtri = Xtri(ind) ; Ytri = Ytri(ind) ;
        % Connectivity
            tri2nod(Elems(tri,:),tri) = 1 ;
        % Shape functions
            % Moments
                Areas(tri) = 1/2*det([ones(3,1) Xtri Ytri]) ;
                aa = cross(Xtri,Ytri) ;
                bb = circshift(Ytri,2) - circshift(Ytri,1) ;
                cc = circshift(Xtri,1) - circshift(Xtri,2) ;
        % Deformation gradient
            % Element-wise
                iiiD(tri,:) = tri ;
                jjjD(tri,:) = Elems(tri,:) ;
                vvvD1(tri,:) = 1/2/Areas(tri)*bb ;
                vvvD2(tri,:) = 1/2/Areas(tri)*cc ;
                %gradT(tri,Elems(tri,:),1,1) = 1/2/Areas(tri)*bb ; %dUx_dx
                %gradT(tri,Elems(tri,:),1,2) = 1/2/Areas(tri)*cc ; %dUx_dy
                %gradT(tri,nNodes+Elems(tri,:),2,1) = 1/2/Areas(tri)*bb ; %dUy_dx
                %gradT(tri,nNodes+Elems(tri,:),2,2) = 1/2/Areas(tri)*cc ; %dUy_dy
        % waitbar
            if toc(ticTime)>0.02
                wtbr = waitbar(tri/nElems,wtbr) ;
                ticTime = tic ;
            end
    end
    wtbr = waitbar(0,wtbr,'Sparsing...') ;
    % TRANSFORM SOME OBJECTS TO SPARSE MATRICES
        D1 = sparse(iiiD(:),jjjD(:),vvvD1(:),nElems,nNodes) ;
        D2 = sparse(iiiD(:),jjjD(:),vvvD2(:),nElems,nNodes) ;
        O = D1*0 ;
        gradT = {[D1 O] [D2 O] ; [O D1] [O D2]} ;
        %gradT = {sparse(gradT(:,:,1,1)) sparse(gradT(:,:,1,2)) ; sparse(gradT(:,:,2,1)) sparse(gradT(:,:,2,2))} ;
        tri2nod = sparse(tri2nod) ;
        Valance = sum(tri2nod,2) ;
        invValance = sparse(1:nNodes,1:nNodes,1./Valance,nNodes,nNodes) ; % Inverse valance (for elemt averaging)
    
    delete(wtbr)
    drawnow ;