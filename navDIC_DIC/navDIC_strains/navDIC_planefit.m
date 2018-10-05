function [obj,hd] = navDIC_planefit(obj,hd)

disp('planefit')

% Params
    %d = 15; % il faudra entrer la valeur servant la construction du maillage
    %width_pln_strains = 300; 3*d ;
    num_pts_fit = 20 ; % Number of pts used to fit the plane

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
                        bx = U(pts_fit,1) ;
                        by = U(pts_fit,2) ;
                        Px = A\bx ;
                        Py = A\by ;
                % Strains
                    E(p,1) = Px(1) ;
                    E(p,2) = Py(2) ;
                    E(p,3) = .5*(Py(1)+Px(2)) ;
            end
        % Save Strains
        obj.Strains(:,:,frame) = E(:,:) ;
    end
    
    
    
% Display Object
    disp(obj)
    
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