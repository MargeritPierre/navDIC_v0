function [obj,hd] = navDIC_planefit(obj,hd)

disp('planefit')

% Params
    d = 40; % il faudra entrer la valeur servant la construction du maillage
    width_pln_strains = 8*d ;

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
    
% If it is the first frame
    if frame == 1
        obj.Strains = zeros(size(Pts,1),3) ;
    end
    
% ELSE, Compute Strains Calculation using displacement ether eulerian (relative displacement) or
% Lagrangian (absolute displacement) 

    
disp(['Computation of Strain ' num2str(disp_Mode)]);
    if frame-refFrame > 1
        % Compute Strains
        % Init
            E(1:nPts,:) = NaN ;
            for p = 1:nPts
                    if isnan(U(p,1)) ; continue ; end ;
                % Get neightoring Pts
                    pt = PtsMov(p,:) ;
                    r = sqrt(sum((PtsMov-repmat(pt,[nPts 1])).^2,2)) ;
                    indFit = r<width_pln_strains ;
                    nfit = size(indFit(indFit),1) ;
                    if nfit<3 ; continue ; end ;
                % Fit a plane
                    %fitPts = Pts(indFit,:) ;
                    xy = PtsMov(indFit,:)-repmat(pt,[nfit 1]) ;
                    % Plane : a+bx+cy = Ux ; d+ex+fy = Uy ; 
                        A = [xy(:,1) xy(:,2) ones(nfit,1)] ;
                        bx = U(indFit,1) ;
                        by = U(indFit,2) ;
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