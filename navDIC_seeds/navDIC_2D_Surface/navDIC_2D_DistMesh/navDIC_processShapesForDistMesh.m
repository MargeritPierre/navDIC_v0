function H = navDIC_processShapesForDistMesh(H)

    % Parameters
        nPtsInt = 10000 ; % num of pts for the intersection computation
        minH0 = 1 ; % default mean bar length
        maxH0 = 200 ;
        minFixCornerAngle = 45*pi/180 ; % minimum angle between two polygon edges to fix a point

    % Is there geometries to process ?
        mesh = [] ;
        if isempty(H.Geometries) ...
                || ~any([H.Geometries.isValid]) ...
            triMesh = findobj(H.Axes,'tag','DistMeshPreview') ;
            delete(triMesh) ;
            return ; 
        end
        
    % Retrieve custom Data
        customData = H.uiContextMenuData(1:length(H.Geometries)) ;
        % Convert to double
            for s = 1:length(customData)
                customData(s).h0 = str2num(customData(s).h0) ;
                customData(s).h = str2num(customData(s).h) ;
                customData(s).l = str2num(customData(s).l) ;
            end
        % If no particular shape has called, reset to current value of h0
            if H.numShape<=0
                h0 = customData(1).h0 ;
            else
                h0 = customData(H.numShape).h0 ;
            end
            for s = 1:length(customData)
                H.uiContextMenuData(s) = setfield(H.uiContextMenuData(s),'h0',num2str(h0)) ;
            end
        % Set variable bar length
            h = [customData.h] ;
            d = [customData.l] ;
            
    % Recursive Distance function
        is1DMesh = true ; % is the mesh 1D curvilinear ?
        dF = {} ; % Individual shape function
        hF = {} ; % Individual shape function
        dFrecurs = {} ; % Recursive dist function
        hFrecurs = {@(p)huniform(p)*h0} ; % Recursive bar length function
        pFix = [] ; % Fixed points
        pInt = [] ; % Intersection Points
        edgPt = {} ;
        geo = 0 ;
        for s = 1:length(H.Shapes)
            if ~H.Geometries(s).isValid ; continue ; end
            geo = geo+1 ;
            pos = H.Geometries(s).Position ;
            % Geometry-dependent properties
                switch H.Geometries(s).Class
                    case 'imrect'
                        is1DMesh = false ;
                        % Local Distance function
                            dF{geo} = @(p)drectangle(p,pos(1),pos(1)+pos(3),pos(2),pos(2)+pos(4)) ;
                        % Trivial Fixed Points
                            pts = [ pos(1) , pos(2) ;...
                                    pos(1) pos(2)+pos(4) ;...
                                    pos(1)+pos(3) pos(2)+pos(4) ;...
                                    pos(1)+pos(3) pos(2) ;...
                                    ] ;
                            pFix = [pFix ; ...
                                    pts ;...
                                    ] ;
                        % Parametrization
                            l = sqrt(sum(diff(pts([1:end,1],:),1,1).^2,2)) ;
                            L = cumsum(l) ;
                            edgPt{s} = @(t)interp1([0;L]/L(end),pts([1:end,1],:),t) ;
                    case 'imellipse'
                        is1DMesh = false ;
                        % Ellipse properties
                            a = pos(3)/2 ;
                            b = pos(4)/2 ;
                            cx = pos(1)+a ;
                            cy = pos(2)+b ;
                        % Distance Function
                            dF{geo} = @(p)(((p(:,1)-cx)./a).^2+((p(:,2)-cy)./b).^2-1)*sqrt(a^2+b^2) ;
                        % Parametrization
                            edgPt{s} = @(t)[cx+a*cos(2*pi*t'),cy+b*sin(2*pi*t')] ;
                    case 'impoly'
                        if ~strcmp(H.Geometries(s).Bool,'impolyline')
                            is1DMesh = false ;
                            % It is a closed polygon
                            polyPos = pos([1:end,1],:) ;
                        else 
                            polyPos = pos ;
                        end
                        % Distance Function
                            dF{geo} = @(p)dpoly(p,polyPos) ;
                        % Fixed Points
                            % Edges vectors
                                edges = diff(polyPos,1) ;
                                vecs = edges./sqrt(sum(edges.^2,2)) ;
                            % Angles at corners
                                angles = acos(sum(vecs.*circshift(vecs,1,1),2)) ;
                            % Fixed Points
                                if ~strcmp(H.Geometries(s).Bool,'impolyline')
                                    % It is a closed polygon
                                    isFixed = angles>minFixCornerAngle ;
                                else 
                                    isFixed = [true ; angles(2:end)>minFixCornerAngle ; true] ;
                                end
                                pFix = [pFix ; pos(isFixed,:)] ;
                        % Parametrization
                            l = sqrt(sum(diff(polyPos,1,1).^2,2)) ;
                            L = cumsum(l) ;
                            edgPt{s} = @(t)interp1([0;L]/L(end),polyPos,t) ;
                    case 'impoint'
                            pFix = [pFix ; pos] ;
                            dF{geo} = @(p)sqrt((p(:,1)-pos(1)).^2+(p(:,2)-pos(2)).^2) ;
                            edgPt{s} = @(t)pos+t(:)*0 ;
                    otherwise
                        disp(H.Geometries(s).Class)
                end
            % Add or Remove Domain
                switch H.Geometries(s).Bool
                    case '-'
                        if geo==1 
                            dFrecurs{geo} = @(p)dF{geo}(p) ;
                        else
                            dFrecurs{geo} = @(p)ddiff(dFrecurs{geo-1}(p),dF{geo}(p)) ;
                        end
                    otherwise %case '+'
                        if geo==1 
                            dFrecurs{geo} = @(p)dF{geo}(p) ;
                        else
                            dFrecurs{geo} = @(p)dunion(dFrecurs{geo-1}(p),dF{geo}(p)) ;
                        end
                end
            % Density function
                hF{geo} = @(p)max(0,1-abs(dF{geo}(p)/d(s))) ;
                hFrecurs{geo+1} = @(p)hFrecurs{geo}(p).*(1-hF{geo}(p))+h(s).*hF{geo}(p) ;
            % Add Fixed Points at Intersections of edges
                if geo>1 && ismember(H.Geometries(s).Class,{'imrect','impoly','imellipse'})
                    t0 = (0:nPtsInt-1)/nPtsInt ;
                    p0 = edgPt{s}(t0) ;
                    d0 = dFrecurs{geo-1}(p0) ;
                    crossEdg = sign(d0)~=circshift(sign(d0),1,1) ;
                    indPt = find(crossEdg) ;
                    for p = 1:length(indPt)
                        tInt = fzero(@(t)dFrecurs{geo-1}(edgPt{s}(t)),t0(indPt(p))) ;
                        pInt(end+1,:) = edgPt{s}(tInt) ;
                    end
                end
        end
    % Final distance Function
        fd = dFrecurs{end} ;
    % Final bar length Function
        FH = @(p)hFrecurs{end}(p) ;
        fh = @(p)FH(p) ; %min(max(FH(p),huniform(p)*minH0),maxH0*huniform(p)) ;
    % Add intersection points to fixed ppoints
        pFix = [pFix ; pInt] ; [] ;
        
    % Remove fixed points that are not on the boundary
        if ~isempty(pFix)
            pFix = pFix(abs(fd(pFix))<=eps,:) ;
        end
        
    % BoundingBox
        margin = 1 ;
        [j,i] = find(H.ROI) ;
        if isempty(j) ; i = pFix(:,1) ; j = pFix(:,2) ; end
        bboxROI = [min(i)-margin min(j)-margin ; max(i)+margin max(j)+margin] ;
        
    % Bar Length Distribution
        if 0 % Density map for debug
            [jj,ii] = meshgrid(1:size(H.ROI,2),1:size(H.ROI,1)) ;
            pij = [jj(:),ii(:)] ;
            hh = fh(pij) ;
            hh(fd(pij)>0) = NaN ;
            hh = reshape(hh,size(ii)) ;
            img = imagesc(repmat(hh,[1 1 1]),'tag','DistMeshPreview') ;
            img.AlphaData = double(~isnan(hh)) ;
            caxis
            pause(1) ;
            delete(findobj(H.Figure,'tag','DistMeshPreview'))
            return ;
        end
        
    % Initial density of points
        D0 = min([h(:);h0]) ;
        
    % Compute the mesh
        if is1DMesh
            mesh = navDIC_computeDistMesh1D(edgPt,fh,pFix,H.Axes) ;
        else
            mesh = navDIC_computeDistMesh2D(fd,fh,D0,bboxROI,pFix,H.Axes) ;
        end
        
    % Return the new Handle structure
        H.DistMesh = mesh ;
        
end
        