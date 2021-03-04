%% STRAIN CRITERION
    
% Gradient matrix, [dUx_dx,dUy_dx,dUx_dy,dUy_dy]
    G = sparse(cat(1,gradT{1,1},gradT{2,1},gradT{1,2},gradT{2,2})) ;

% Linearized Green-Lagrane Epsilon = BB*a, epsilon 
    B = sparse(cat(1,gradT{1,1},gradT{2,2},0.5*(gradT{1,2}+gradT{2,1}))) ;
    BB = B'*B ; % To minimize the "strain energy"
    
% Minimisation of the strain GRADIENT in the dir. of the edge normal
    % Build edges associated to elements
        elmtEDGES = [Elems(:,1:2);Elems(:,2:3);Elems(:,[3,1])] ;
        elmtEDGES = sort(elmtEDGES,2) ;
        elmtEDGES = elmtEDGES(:,1)+1i*elmtEDGES(:,2) ; % Complex notation is more convenient
        elmtEDGES = reshape(elmtEDGES,[nElems 3]) ;
    % Edges of the mesh
        EDGES = unique(elmtEDGES) ;
        nEdges = length(EDGES) ;
    % Edge Connectivity
        edg2nod = sparse(nNodes,nEdges) ; % nodes linked to each edge
        tri2edg = sparse(nEdges,nElems) ; % elements linked to each edge
    % Edge constraint
        switch strainCriterion
            case 'normal'
                Ed = sparse(2*nEdges,4*nElems) ; % projection on the normal
            case 'full'
                Ed = sparse(4*nEdges,4*nElems) ; % full gradient
        end
    % Process
        wtbr = waitbar(0,'Edge constraint...') ;
        for edg = 1:nEdges
            % Retrieve
                thisEdge = [real(EDGES(edg)) imag(EDGES(edg))] ;
                edg2nod(thisEdge,edg) = 1 ;
            % Find the attached triangle(s)
                [elmts,~] = find(abs(elmtEDGES-EDGES(edg))<eps) ;
                tri2edg(edg,elmts) = 1 ;
                if length(elmts)==1 ; continue ; end % the edge is naked
            % Edge normal
                tang = diff(Nodes(thisEdge,:),1,1) ; tang = tang.'/norm(tang) ;
                normal = [tang(1) tang(2) ; -tang(2) tang(1)]\[0;1] ;
                dist = diff([Xc(elmts(:))';Yc(elmts(:))'],1,2) ;
                dVect = [-1 1]/norm(dist) ;
            % Constrain
                switch strainCriterion
                    case 'normal' % Projection of the STRAINS on the normal
                        % dExx.nx + dExy.ny = 0
                            Ed(edg,elmts) = normal(1)*dVect ;
                            Ed(edg,elmts+2*nElems) = normal(2)*dVect ;
                        % dExy.nx + dEyy.ny = 0
                            Ed(edg+nEdges,elmts+2*nElems) = normal(1)*dVect ;
                            Ed(edg+nEdges,elmts+1*nElems) = normal(2)*dVect ;
                    case 'full' % variation of the GRADIENT between elements
                        for comp = 1:4
                            Ed(edg+(comp-1)*nEdges,elmts+(comp-1)*nElems) = dVect ;
                        end
                end
            % Waitbar
                wtbr = waitbar(edg/nEdges,wtbr) ;
        end
        delete(wtbr) ;
    % Sparsify...
        Ed = sparse(Ed) ; % gives the difference between two elements
        Em = double(logical(abs(Ed))) ;
        edg2nod = sparse(edg2nod) ;
        tri2edg = sparse(tri2edg) ;