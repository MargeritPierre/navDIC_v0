function mesh = navDIC_computeDistMesh1D(edgFcn,fh,pfix,Axes)
%DISTMESH1D 1-D Mesh Generator 
    
% Delete previous graphic handles
    meshTag = 'DistMeshPreview' ;
    gHandles = findobj(Axes,'tag',meshTag) ;
    delete(gHandles) ;

% Build the mesh
    p = [] ; bars = [] ; t0 = linspace(0,1,1000) ;
    for s = 1:numel(edgFcn)
    % initial uniform distribution
        p0 = edgFcn{s}(t0) ; % initial points
        dl = sqrt(sum(diff(p0,1,1).^2,2)) ; % initial spacing
        dl0 = fh(.5*(p0(1:end-1,:)+p0(2:end,:))) ; % local target spacing
    % number of points
        n = dl./dl0 ; % number of points in each interval
        n = cumsum([0 ; n]) ; % cummulative number
        n = (n./n(end)).*floor(n(end)) ; % integer number
    % new line parameter
        if n(end)>0
            t = interp1(n,t0,(0:n(end))') ; % parameters
        else
            t = [0;1] ;
        end
    % new points
        ep = edgFcn{s}(t(:)) ;
        [uep,ia,ic] = unique(ep,'rows') ;
    % new bars
        ebars = size(p,1) + [1:numel(t)-1 ; 2:numel(t)]' ;
        ebars = reshape(ic(ebars),size(ebars)) ;
        ebars(ebars(:,1)==ebars(:,2),2) = NaN ;
    % update
        bars = [bars ; ebars] ;
        p = [p ; uep] ;
    end
    
% Plot
    triMesh = patch('vertices',p...
                    ,'faces',padarray(bars,[0 1],NaN,'post') ...
                    ,'edgecolor','b' ...
                    ,'marker','.' ...
                    ,'markersize',20 ...
                    ,'facecolor','none'...
                    ,'tag',meshTag...
                    ...,'hittest','off' ...
                    ,'pickableparts','none' ...
                    ) ;
    uistack(triMesh,'bottom') ;
    uistack(findobj(Axes,'type','image'),'bottom')
    
% DISPLAY INFOS
    infos = {} ;
    infos{end+1} = {[num2str(size(p,1)),' nodes']} ;
    infos{end+1} = {[num2str(size(t,1)),' triangles']} ;
    disp(newline)
    disp('--- DistMesh 1D ---')
    for in = [infos{:}]
        disp(['  ',in{1}])
    end
    disp('-------------------')
    
% Set the mesh Callbacks
    %setCallbacks() ;
    
% Mesh to return
    mesh.Points = p ;
    mesh.Triangles = bars ;
    mesh.Patches = triMesh ; %patch(H.Axes,'vertices',p,'faces',bars) ;
    
    
% ===================================================================================================================    
% CHILD FUNCTIONS
% ===================================================================================================================

   
    function setCallbacks()
        % Find the element which is under attention
            fig = gcf ;
            ax = gca ;
            pl = plot(ax,NaN,NaN,'color',triMesh.EdgeColor,'hittest','off','tag',meshTag,'pickableparts','none') ;
            mousePos = @()ax.CurrentPoint(1,1:2) ;
            fig.WindowButtonMotionFcn = @(src,evt)showSelectedElements(pl,mousePos(),triMesh) ;
            fig.WindowButtonDownFcn = @(src,evt)clickCallback(mousePos(),triMesh,evt) ;
    end

    function selected = selectedElement(mp,p,t)
        selected = [] ;
        if isempty(t) ; return ; end
        x = p(:,1) ; y = p(:,2) ;
        xc = mean(x(t),2) ; yc = mean(y(t),2) ; 
        distToElems = sqrt((xc-mp(1)).^2 + (yc-mp(2)).^2) ;
        [~,closestElems] = sort(distToElems,'ascend') ;
        for e = 1:min(10,numel(t(:,1)))
            elmt = closestElems(e) ;
            xe = x(t(elmt,:)) ; ye = y(t(elmt,:)) ; 
            if inpolygon(mp(1),mp(2),xe,ye) ; selected(end+1) = elmt ; end
        end
    end

    function showSelectedElements(pl,mp,mesh)
        if ~isvalid(mesh) ; return ; end
        p = mesh.Vertices ; t = mesh.Faces ;
        selected = selectedElement(mp,p,t) ;
        if isempty(selected)
            pl.XData = NaN ; 
            pl.YData = NaN ; 
            return ; 
        end
        selected = selected(1) ;
        x = p(:,1) ; y = p(:,2) ;
        xe = x(t(selected,[1:end,1])) ; ye = y(t(selected,[1:end,1])) ;
        if ~isequal(pl.XData,xe') && ~isequal(pl.YData,ye') 
            pl.XData = xe' ;
            pl.YData = ye' ;
            disp(['Distmesh: element(s) ',mat2str(selected)])
        end
    end

    function clickCallback(mp,mesh,evt)
        switch evt.Source.SelectionType
            case 'normal' % Left-Click
            case 'alt' % Right-Click
                selected = selectedElement(mp,mesh.Vertices,mesh.Faces) ;
                %mesh.Faces(selected,:) = [] ;
            case 'extend' % Scroll Wheel
        end
    end


end
