function [obj,hd] = navDIC_lineStrain(obj,hd)

disp('lineStrain')

% Get the line direction (from end points)
    v0 = obj.Points(end,:)-obj.Points(1,:) ; v0 = v0/norm(v0) ;
    v = obj.MovingPoints(end,:,hd.CurrentFrame) - obj.MovingPoints(1,:,hd.CurrentFrame) ; v = v/norm(v) ;
    
% Positions in the line coordinate
    x0 = obj.Points*v0(:) ; x0 = x0-x0(1) ;
    x = obj.MovingPoints(:,:,hd.CurrentFrame)*v(:) ; % x = x-x(1) ;
    dx = x-x0 ;
    
% Valid indices
    valid = ~isnan(dx) ;
    dx = dx(valid) ;
    x0 = x0(valid) ;
    
% Least-squares fit dx = a * x0 + b
    ab = [x0(:) ones(length(x0),1)]\dx(:) ;
    obj.Strains(:,:,hd.CurrentFrame) = ab(1) ;
    
    
    
% disp('deltaL_L')
%     
% % Retrieve Infos
%     frame = hd.CurrentFrame ;
%     camID = obj.CamIDs ;
%     PtsMov = round(obj.MovingPoints(:,:,frame)) ; % 
%     disp_Mode = obj.displMode ;
%     U = obj.Displacements ;
%     Pts = obj.Points; % ref absolu dl/l0
%     nPts = length(Pts(:,:)) ;
%     L0 = ones(nPts-1,1) ;
%     if strcmpi(disp_Mode,'abs') || frame==1 
%         for i=1:(nPts-1)
%             L0(i) = (Pts(i,:) * Pts(i+1,:)') ^ 0.5 ;
%         end
%     else
%         PtsMovRef = round(obj.MovingPoints(:,:,frame-1)); % ref relatif dl/l 
%         for i=1:(nPts-1)
%             L0(i) = (PtsMovRef(i,:) * PtsMovRef(i+1,:)') ^ 0.5 ;
%         end
%     end
%         
% % If it is the first frame
%     if frame == 1
%         obj.Strains = zeros(nPts-1,1) ;
%     end
%     
% % ELSE, Compute Strains Calculation using displacement ether eulerian (relative displacement) or
% % Lagrangian (absolute displacement) 
% 
% disp(['Computation of Strain ' num2str(disp_Mode)]);
% if frame > 1
% 
%     % Compute Strains
%     % Init
%         L(1:nPts-1,1) = NaN ;
%         E(1:nPts-1,1) = NaN ;
%         for p = 1:(nPts-1)
%             if isnan(U(p+1,1)) || isnan(U(p,1))  ; continue ; end ;
%         % Actual length
%             L(p) = (PtsMov(p+1,:) * PtsMov(p,:)') ^ .5 ;
%         % Strains
%             E(p) = (L(p)-L0(p)) / L0(p) ;
%         end
%     % Save Strains
%     obj.Strains(1,:,frame) = E(:) ;
% end
% 
% % Display Object
%     disp(obj)