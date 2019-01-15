function [obj,hd] = navDIC_planefit(obj,hd)

%disp('planefit')
startTime = tic() ;

% Params
    %d = 15; % il faudra entrer la valeur servant la construction du maillage
    %width_pln_strains = 300; 3*d ;
    num_pts_fit = 100 ; % Number of pts used to fit the plane
    FIT = 'TLS' ; % Type of plane-fit (see below)

% Retrieve Infos
    frame = hd.CurrentFrame ;
    camID = obj.CamIDs ;
    refFrame = obj.RefFrame ;
    U = obj.Displacements(:,:,frame) ;
    Pts = obj.Points;
    try
        PtsMov = round(obj.MovingPoints(:,:,frame));
    catch
        PtsMov = Pts;
    end
    disp_Mode = obj.displMode ;
    nPts = length(PtsMov(:,:));
    
    % Blah blah
    if frame<=obj.RefFrame
        obj.Strains = zeros([size(Pts,1) 3 frame]) ;
    end
    
% ELSE, Compute Strains Calculation using displacement ether eulerian (relative displacement) or
% Lagrangian (absolute displacement) 

    
disp(['Computation of Strain ' num2str(disp_Mode)]);
    if frame-refFrame >= 1
        % Compute Strains
        % Init
            E(1:nPts,1:3) = NaN ;
            for p = 1:nPts
                    if isnan(U(p,1)) ; continue ; end ;
                % Get neightoring Pts
                    pt = PtsMov(p,:) ;
                    r = sqrt(sum((PtsMov-repmat(pt,[nPts 1])).^2,2)) ;
                    [~,indsort] = sort(r) ;
                    pts_fit = indsort(1:num_pts_fit) ;
                % Fit a plane
                    xy = PtsMov(pts_fit,:)-repmat(pt,[num_pts_fit 1]) ;
                    % Plane : a+bx+cy = Ux ; d+ex+fy = Uy ; 
                        A = [xy(:,1) xy(:,2) ones(num_pts_fit,1)] ;
                        b = U(pts_fit,:) ;
                        switch FIT
                            case 'LS' % LEAST SQUARES
                                P = A\b ;
                            case 'TLS' % TOTAL LEAST SQUARES
                                [~,~,v] = svd([A b],0);           % find the SVD of Z.
                                vXY = v(1:3,4:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
                                vYY = v(4:end,4:end) ; % Take the bottom-right block of V.
                                P = -vXY/vYY;
                            case 'GLS' % GENERALIZED LEAST SQUARES
                                P = A\b ;
                                res  = sum((A*P-b).^2,2)/2 ;
                                RES = diag(1./res) ;
                                P = (RES*A)\(RES*b) ;
                        end
                % Strains
                    E(p,1) = P(1,1) ;
                    E(p,2) = P(2,2) ;
                    E(p,3) = .5*(P(1,2)+P(2,1)) ;
            end
        % Save Strains
        obj.Strains(:,:,frame) = E(:,:) ;
        
    end
    
    % TIMING
        disp(['planefit: ',num2str(toc(startTime)*1000,'%.1f'),' ms']) ;
    
    
    
% Display Object
    %disp(obj)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD
%             % Update Plots
%             if exist('pl')
%                 pl.XData = Pts(:,1) ;
%                 pl.YData = Pts(:,2) ;
%             end
%             im.CData = repmat(imgs{i},[1 1 3]) ;
%             srf.Vertices = [Pts zeros(nPts,1)] ;
%             srf.FaceVertexCData = E(:,2,i)*100 ; sqrt(sum(U.^2,2)) ;
%             drawnow ;
%             pause(.01) ;