function H = navDIC_processShapesForDistMesh(H)

    % Parameters
        nPtsInt = 300 ; % num of pts for the intersection computation
        minH0 = 1 ; % default mean bar length
        maxH0 = 200 ;

    % Is there geometries to process ?
        mesh = [] ;
        if isempty(H.Geometries) || ~any([H.Geometries.isValid]) 
            triMesh = findobj(H.Axes,'tag','DistMeshPreview') ;
            delete(triMesh) ;
            return ; 
        end
        
    % Retrieve custom Data
        customData = H.uiContextMenuData ;
        % Convert to double
            for s = 1:length(customData) ;
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
            for s = 1:length(customData) ;
                H.uiContextMenuData(s) = setfield(H.uiContextMenuData(s),'h0',num2str(h0)) ;
            end
        % Set variable bar length
            h = [customData.h] ;
            d = [customData.l] ;
            
    % Recursive Distance function
        dF = {} ; % Individual shape function
        hF = {} ; % Individual shape function
        dFrecurs = {} ; % Recursive dist function
        hFrecurs = {@(p)huniform(p)*h0} ; % Recursive bar length function
        pFix = [] ; % Fixed points
        pInt = [] ; % Intersection Points
        geo = 0 ;
        for s = 1:length(H.Shapes)
            if ~H.Geometries(s).isValid ; continue ; end
            geo = geo+1 ;
            pos = H.Geometries(s).Position ;
            % Geometry-dependent properties
                switch H.Geometries(s).Class
                    case 'imrect'
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
                            edgPt = @(t)interp1([0;L]/L(end),pts([1:end,1],:),t) ;
                    case 'imellipse'
                        % Ellipse properties
                            a = pos(3)/2 ;
                            b = pos(4)/2 ;
                            cx = pos(1)+a ;
                            cy = pos(2)+b ;
                        % Distance Function
                            dF{geo} = @(p)(((p(:,1)-cx)./a).^2+((p(:,2)-cy)./b).^2-1)*sqrt(a^2+b^2) ;
                        % Parametrization
                            edgPt = @(t)[cx+a*cos(2*pi*t'),cy+b*sin(2*pi*t')] ;
                    case 'impoly'
                        % Distance Function
                            dF{geo} = @(p)dpoly(p,pos([1:end,1],:)) ;
                        % Fixed Points
                            pFix = [pFix ; pos] ;
                        % Parametrization
                            l = sqrt(sum(diff(pos([1:end,1],:),1,1).^2,2)) ;
                            L = cumsum(l) ;
                            edgPt = @(t)interp1([0;L]/L(end),pos([1:end,1],:),t) ;
                end
            % Add or Remove Domain
                switch H.Geometries(s).Bool
                    case '+'
                        if geo==1 
                            dFrecurs{geo} = @(p)dF{geo}(p) ;
                        else
                            dFrecurs{geo} = @(p)dunion(dFrecurs{geo-1}(p),dF{geo}(p)) ;
                        end
                    case '-'
                        if geo==1 
                            dFrecurs{geo} = @(p)dF{geo}(p) ;
                        else
                            dFrecurs{geo} = @(p)ddiff(dFrecurs{geo-1}(p),dF{geo}(p)) ;
                        end
                end
            % Density function
                hF{geo} = @(p)max(0,1-abs(dF{geo}(p)/d(s))) ;
                hFrecurs{geo+1} = @(p)hFrecurs{geo}(p).*(1-hF{geo}(p))+h(s).*hF{geo}(p) ;
            % Add Fixed Points at Intersections of edges
                if geo>1
                    t0 = (0:nPtsInt-1)/nPtsInt ;
                    p0 = edgPt(t0) ;
                    d0 = dFrecurs{geo-1}(p0) ;
                    crossEdg = sign(d0)~=circshift(sign(d0),1,1) ;
                    indPt = find(crossEdg) ;
                    for p = 1:length(indPt)
                        tInt = fzero(@(t)dFrecurs{geo-1}(edgPt(t)),t0(indPt(p))) ;
                        pInt(end+1,:) = edgPt(tInt) ;
                    end
                end
        end
    % Final dist. Function
        fd = dFrecurs{end} ;
    % Final dist. Function
        FH = @(p)hFrecurs{end}(p) ;
        fh = @(p)FH(p) ; %min(max(FH(p),huniform(p)*minH0),maxH0*huniform(p)) ;
    % Add intersection points to fixed ppoints
        pFix = [pFix ; pInt] ;
        
    % Remove fixed points that are not on the boundary
        if ~isempty(pFix)
            pFix = pFix(abs(fd(pFix))<=eps,:) ;
        end
        
    % BoundingBox
        [j,i] = find(H.ROI) ;
        margin = 1 ;
        bboxROI = [min(i)-margin min(j)-margin ; max(i)+margin max(j)+margin] ;
        
    % Bar Length Distribution
        [jj,ii] = meshgrid(1:size(H.ROI,2),1:size(H.ROI,1)) ;
        pij = [jj(:),ii(:)] ;
        hh = fh(pij) ;
        hh(fd(pij)>0) = NaN ;
        if 0 % Density map for debug
            hh = reshape(hh,size(ii)) ;
            img = imagesc(repmat(hh,[1 1 1]),'tag','DistMeshPreview') ;
            img.AlphaData = double(~isnan(hh)) ;
            caxis
            pause(1) ;
            delete(findobj(H.Figure,'tag','DistMeshPreview'))
            return ;
        end
    % Initial density of points 
        D0 = min(min(hh(~isnan(hh)))) ;
        
    % Compute the mesh
        mesh = navDIC_computeDistMesh2D(fd,fh,D0,bboxROI,pFix,H.Axes) ;
        
    % Return the new Handle structure
        H.DistMesh = mesh ;
        
end
        