% Integral over a circular domain

global hd

crackTipPosition = hd.Seeds(3).MovingPoints ;
energySeed = hd.Seeds(2) ;
IMG = repmat(hd.Images{1},[1 1 3 1]) ;
dataField = 'W' ;
domainRadius = 50 ;
domainShift = [0 -1]*domainRadius ;
nPint = 2500 ;
config = 'reference' ; 'current' ;
showZone = true ;

% Create a domain in which to evaluate the integral
    the = linspace(0,2*pi,100)' ; 
    pDom = domainRadius*[cos(the) sin(the)]+domainShift ;

% Discretize the domain (quadrature points)
    margin = 0.1 ;
    domArea = polyarea(pDom(:,1),pDom(:,2)) ;
    domLims = [min(pDom,[],1) ; max(pDom,[],1)]+range(pDom,1).*[-1;1]*margin ;
    domRangeArea = prod(range(domLims,1)) ;
    nPint = round(nPint*domRangeArea/domArea) ;
    da = domRangeArea/nPint ;
    dl = sqrt(da) ;
    xx = domLims(1,1):dl:domLims(2,1) ;
    yy = domLims(1,2):dl:domLims(2,2) ;
    [XX,YY] = meshgrid(xx,yy) ;
    MASK = inpolygon(XX(:),YY(:),pDom(:,1),pDom(:,2)) ;
    XX = XX(MASK) ; YY = YY(MASK) ;
    
% Init the figure
    if showZone
        clf reset ;
        im = imagesc(IMG(:,:,:,1)) ;
            %colormap gray 
            axis ij
            axis equal
            axis tight
            set(gca,'xtick',[],'ytick',[])
        pa = patch('Faces',energySeed.Triangles,'Vertices',[energySeed.MovingPoints(:,:,1) energySeed.MovingPoints(:,1,1)*0],'edgecolor','none','facecolor','interp') ;
        pl = patch(pDom(:,1),pDom(:,2),pDom(:,2)*0,'facecolor','none','edgecolor','k','linewidth',2) ;
    end
    
% Evaluate the integral ;
    INT = NaN(hd.nFrames,1) ;
    for fr = 1:hd.nFrames
        % Crack tip
            cen = crackTipPosition(:,:,fr) ;
            if any(isnan(cen)) ; continue ; end
        % Evaluation points
            switch config
                case 'current'
                    pts = energySeed.MovingPoints(:,:,fr) ;
                    pInt = [XX(:) YY(:)] + cen ;
                case 'reference'
                    pts = energySeed.MovingPoints(:,:,1) ;
                    cen = full(energySeed.interpMat(cen,'extrap',energySeed.MovingPoints(:,:,fr)))*pts ;
                    pInt = [XX(:) YY(:)] + cen ;
            end
        % Evaluate
            T = energySeed.interpMat(pInt,0,pts) ;
            field = energySeed.DataFields.(dataField)(:,fr) ;
            data = T*field ;
            INT(fr) = sum(data)*da/domArea ;
        % Display
            if showZone
                switch config
                    case 'current'
                        im.CData = IMG(:,:,:,fr) ;
                        pa.Vertices(:,1:2) = pts ;
                    case 'reference'
                end
                pa.FaceVertexCData = field ;
                pl.Vertices = pDom + cen ;
                drawnow ;
            end
    end

% Plot the result
     if showZone ; clf ; end
        plot(INT)
        ax = gca ;
        axis auto
        ax.YLim(1) = 0 ;




