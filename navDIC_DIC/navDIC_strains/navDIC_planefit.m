function [obj,hd] = navDIC_planefit(obj,hd)

%disp('planefit')

% Params
    %d = 15; % il faudra entrer la valeur servant la construction du maillage
    %width_pln_strains = 300; 3*d ;
    num_pts_fit = 30 ; % Number of pts used to fit the plane (can be INF)
    FIT = 'TLS' ; % Type of plane-fit (see below)
    expr = 'full' ; % Green strains: 'linearized' or 'full'

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    refFrame = obj.RefFrame ;
    disp_Mode = obj.displMode ;
    U = obj.Displacements(:,:,frame) ;
    Pts = obj.Points;
    PtsMov = Pts ; obj.MovingPoints(:,:,frame);
    PtsMov(isnan(obj.MovingPoints(:,:,frame))) = NaN ;
    nPts = length(PtsMov(:,:));
    
% Modify parameters
    if num_pts_fit>nPts ; num_pts_fit = nPts ; end
    
% Before the reference frame
    if frame<=obj.RefFrame
        obj.Strains = zeros([size(Pts,1) 3 frame]) ;
        return
    end
    
% Compute Strains 
    % Init
        E(1:nPts,1:3) = NaN ;
        startTime = tic() ;
    % Loop over the points
        for p = 1:nPts
            % If the point is not valid, skip
                if any(isnan(U(p,:))) ; continue ; end
            % Get neightoring Pts
                pt = PtsMov(p,:) ;
                r = sqrt(sum((PtsMov-repmat(pt,[nPts 1])).^2,2)) ;
                [~,indsort] = sort(r) ;
                pts_fit = indsort(1:num_pts_fit) ;
            % Cull NaNs
                pts_fit = pts_fit(~any(isnan(PtsMov(pts_fit,:)),2)) ;
            % Fit a plane
                % Get positions
                    xy = bsxfun(@minus,PtsMov(pts_fit,:),pt) ;
                % Plane : a+bx+cy = Ux ; d+ex+fy = Uy ; 
                    A = [xy(:,1) xy(:,2) ones(length(pts_fit),1)] ;
                    b = U(pts_fit,:) ;
                    switch FIT
                        case 'LS' % LEAST SQUARES
                            P = A\b ;
                        case 'TLS' % TOTAL LEAST SQUARES
                            [~,~,v] = svd([A b],0) ;           % find the SVD of Z.
                            vXY = v(1:3,4:end) ;     % Take the block of V consisting of the first n rows and the n+1 to last column
                            vYY = v(4:end,4:end) ; % Take the bottom-right block of V.
                            P = -vXY/vYY ;
                        case 'GLS' % GENERALIZED LEAST SQUARES
                            P = A\b ;
                            res  = sum((A*P-b).^2,2)/2 ;
                            RES = diag(1./res) ;
                            P = (RES*A)\(RES*b) ;
                    end
            % Strains
                switch expr
                    case 'linearized'
                        E(p,1) = P(1,1) ;
                        E(p,2) = P(2,2) ;
                        E(p,3) = .5*(P(1,2) + P(2,1)) ;
                    case 'full'
                        E(p,1) = P(1,1) + 0.5*(P(1,1)^2+P(1,2)^2) ;
                        E(p,2) = P(2,2) + 0.5*(P(2,2)^2+P(2,1)^2) ;
                        E(p,3) = .5*(P(1,2) + P(2,1) + P(1,1)*P(2,1) + P(2,2)*P(1,2)) ;
                end
            % If the number of neightbors is equal to the number of points
                if num_pts_fit == nPts
                    E = repmat(E(p,:),[nPts 1]) ;
                    break ;
                end
        end
    % Save Strains
        obj.Strains(:,:,frame) = E ;
    % Timing
        disp(['planefit (',num2str(disp_Mode),'):',num2str(toc(startTime)*1000,'%.1f'),' ms']) ;
    