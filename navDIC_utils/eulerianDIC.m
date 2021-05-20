%% GLOBAL DIC IN EULERIAN FRAME
% Before that, please create a gridMeshSeed
%open GridMeshSeed

%% LOAD AND PREPARE THE SEED

    global hd
    seed = hd.Seeds(1) ;
    
    % External objects
    mesh = pkg.geometry.mesh.Mesh('Nodes',seed.Points,'Elems',seed.Elems) ;
    refImg = seed.refImgs{end} ;
    
    % Region of Interest
    [N,INSIDE] = pkg.dic.global.mapping(refImg,mesh) ;
    ROI = find(any(INSIDE,2)) ;
    N = N(ROI,:) ;
    INSIDE = INSIDE(ROI,:) ;
    
    
%% PERFORM DIC
    switch 3 % Transfer matrix T: constrain DOFs s.t V = T*v, v constrained
        case 1 % All DOFS
            T = speye(2*mesh.nNodes) ;
        case 2 % Uniform X speed by element
            e2n = mesh.elem2node('mean') ;
            T = [e2n ; e2n*0] ;
        case 3 % Uniform X speed by vertical edge
            edgN = mesh.getNormals(mesh.Edges) ; % edge normals
            vEdg = abs(sum(edgN.*[1 0 0],2))>0.9 ; % vertical edges
            e2n = mesh.edge2node ; % edge to node connectivity
            ve2n = e2n(:,vEdg) ; % vertical edge connectivity
            T = [ve2n ; ve2n*0] ; 
            if 0 % add vertical contraction: difference in vertical edge's nodes V2
                el = pkg.data.sparse2list(ve2n') ; % list of nodes attached to each vertical edge
                ye = reshape(mesh.Nodes(el,2),size(el)) ; % node elevation
                sn = round(2*(ye-mean(ye,2))./diff(ye,1,2)) ; % node sign [-1 or 1]
                Mc = sparse(el(:)',repmat(1:size(el,1),[1 2]),sn(:)',mesh.nNodes,size(el,1)) ;
                T = [T [Mc*0 ; Mc]] ;
            end
    end
    filtSize = 3*[1 1] ; % image filter size
    usePreviousVelocity = true ; % init using previously obtained velocity
    useAcceleration = false ; % init using first order acc. extrapolation
    tolV = 1e-2 ; % maximum velocity change at convergence
    maxIt = 100 ; % maximum number of iterations
    approxGradients = true ; % approximate dg_dx with dG_dX
    plotFreq = 0 ; % plot update frequency
    frames = 2:hd.nFrames ;
    
    clf
    im = image(refImg) ;
    plot(mesh,'VisibleFaces','none') ;
    qui = quiver(mesh.Nodes(:,1),mesh.Nodes(:,2),mesh.Nodes(:,1)*0,mesh.Nodes(:,2)*0,'AutoScale','off') ;
    axis ij
    axis equal
    axis tight
    
    X = seed.Points ; 
    V = zeros(mesh.nNodes,2,hd.nFrames) ;
    
    filter = blackman(filtSize(2)+2)'.*blackman(filtSize(1)+2) ;
    filter = filter(2:end-1,2:end-1) ;
    filter = filter/sum(filter(:)) ;
    
    
    lastPlotTime = uint64(0) ;
    for fr = frames
    % Get Images
        G = double(hd.Images{seed.CamIDs}{fr-1}) ;
        g = double(hd.Images{seed.CamIDs}{fr}) ;
    % Filter images
        if any(filtSize>1) 
            G = conv2(G,filter,'same') ;
            g = conv2(g,filter,'same') ;
        end
    % Display images
        im.CData = imfuse(G,g,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    % Current Image Gradients
        if approxGradients
            dG_dX = conv2(G,[1 0 -1]/2,'same') ;
            dG_dY = conv2(G,[1 0 -1]'/2,'same') ;
            dr_dv = [dG_dX(ROI).*N dG_dY(ROI).*N] ; 
            H = dr_dv'*dr_dv ;
        else
            dg_dx = conv2(g,[1 0 -1]/2,'same') ;
            dg_dy = conv2(g,[1 0 -1]'/2,'same') ;
        end
    % Init Velocity
        if usePreviousVelocity 
            valid = all(~isnan(V(:,:,fr-1)),2) ;
            V(valid,:,fr) = V(valid,:,fr-1) ; 
            if useAcceleration && fr>=3
                valid = valid & all(~isnan(V(:,:,fr-2)),2) ;
                V(valid,:,fr) = V(valid,:,fr) + diff(V(valid,:,fr+(-2:-1)),1,3) ;
            end
        end
    % DIC Loop
        outFlag = [] ; it = 0 ; valid = true(mesh.nNodes,1) ;
        while ~any(outFlag)
        % Current image & gradients
            xx = N*(X+V(:,:,fr)) ;
            gi = interp2(g,xx(:,1),xx(:,2)) ;
            if ~approxGradients
                dgi_dx = interp2(dg_dx,xx(:,1),xx(:,2)) ;
                dgi_dy = interp2(dg_dy,xx(:,1),xx(:,2)) ;
            end
        % Residual, Jacobian & hessian
            r = gi-G(ROI) ;
            if ~approxGradients
                dr_dv = [dgi_dx.*N dgi_dy.*N] ; 
                H = dr_dv'*dr_dv ;
            end
            J = dr_dv'*r ;
        % Valid geometry & constaints
            vv = [valid;valid] ;
            Tv = T(vv,:) ;
            Jv = Tv'*J(vv) ;
            Hv = Tv'*H(vv,vv)*Tv ;
        % Update
            dv = Hv\Jv ;
            V(:,:,fr) = V(:,:,fr) - reshape(T*dv,[],2) ;
        % Infos
            disp([ ...
                'Frame ' num2str(fr) ...
                ' it ' num2str(it) ...
                ' max(dV) ' num2str(max(abs(dv(~isnan(dv))))) ...
                ]) ;
        % Break ?
            it = it+1 ;
            valid = all(~isnan(V(:,:,fr)),2) ;
            if it>=maxIt ; outFlag = 'maxIt' ; end
            if max(abs(dv(~isnan(dv))))<tolV ; outFlag = 'tolV' ; end
            if ~any(valid) ; outFlag = 'allNaN' ; end
        % Display
            if any(outFlag) || toc(lastPlotTime)>1/plotFreq
                qui.UData = V(:,1,fr) ;
                qui.VData = V(:,2,fr) ;
                drawnow ;
                lastPlotTime = tic ;
            end
        end
        disp(['outFlag: ' outFlag]) ; 
    end
    %return ;
    
% Compute Data Fields
    DATA = [] ;
    DATA.Position = 'Position' ;
        DATA.NaN = NaN(mesh.nNodes,1,hd.nFrames) ;
        DATA.x1 = repmat(mesh.Nodes(:,1),[1 1 hd.nFrames]) ;
        DATA.x2 = repmat(mesh.Nodes(:,2),[1 1 hd.nFrames]) ;
    DATA.Velocity = 'Velocity' ;
        DATA.v1 = V(:,1,:) ;
        DATA.v2 = V(:,2,:) ;
        DATA.V = sqrt(sum(V.^2,2)) ;
    DATA.StrainRate = 'Strain Rate' ;
        G = mesh.gradMat ;
        DATA.D11 = reshape(G{1}*reshape(V(:,1,:),[],hd.nFrames),[],1,hd.nFrames) ;
        DATA.D22 = reshape(G{2}*reshape(V(:,2,:),[],hd.nFrames),[],1,hd.nFrames) ;
        DATA.D12 = 0.5*reshape(G{1}*reshape(V(:,2,:),[],hd.nFrames)+G{2}*reshape(V(:,1,:),[],hd.nFrames),[],1,hd.nFrames) ;
    DATA.MassFlow = 'Mass Flow' ;
        massThreshold = 43 ; % use imageThreshold to determine this (tSlider.Value)
        img = cat(4,hd.Images{seed.CamIDs}{frames}) ;
        density = double(reshape(img,[],numel(frames))>=massThreshold) ;
        DATA.M = DATA.NaN(1:mesh.nElems,:,:) ; 
        DATA.M(:,frames) = INSIDE'*density(ROI,:) ;
        DATA.Q = DATA.NaN ;
        DATA.Q(:,frames) = DATA.V(:,frames).*(mesh.elem2node('mean')*DATA.M(:,frames)) ;
    
% PUSH TO NAVDIC
    seed.MovingPoints = repmat(seed.Points,[1 1 hd.nFrames]) ;
    seed.computeDataFields ;
    seed.DataFields = DATA ;
    for prev = hd.Previews(:)'
    	prev.updatePreview(hd) ;
    end


