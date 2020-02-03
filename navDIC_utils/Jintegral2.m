%% J-INTEGRAL ESTIMATION ON DIC DATA


%% IMPORT THE SEED AND COMPUTE STRESSES AND STRAIN ENERGY
    seedNumber = 2 ;
    
    % Import hd struct
        global hd
        seed = hd.Seeds(seedNumber) ;
    % Re-compute the seed Data Fields on points
        seed.DataOnNodes = true ;
    % Extract the total strain
        EPS = seed.DataFields.TSeq ;
        %EPS = cumsum(seed.DataFields.Deq,3) ;
        
%% COMPUTE EQUIVALENT CAUCHY STRESSES
    %% Ramberg_Osgood
        [SIG,Es,Nus] = Ramberg_Osgood(EPS) ;
            seed.DataFields.CauchyStresses = 'Cauchy Stresses' ;
                seed.DataFields.S11 = Es./(1-Nus.^2).*(seed.DataFields.TS11 + Nus.*seed.DataFields.TS22) ;
                seed.DataFields.S22 = Es./(1-Nus.^2).*(seed.DataFields.TS22 + Nus.*seed.DataFields.TS11) ;
                seed.DataFields.S12 = Es./(1+Nus).*seed.DataFields.TS12 ;
                seed.DataFields.Seq = SIG ;
                
    %% Power hardening
        [SIG,Es,Nus] = Power_Hardening(EPS) ;
            seed.DataFields.CauchyStresses = 'Cauchy Stresses' ;
                seed.DataFields.S11 = Es./(1-Nus.^2).*(seed.DataFields.TS11 + Nus.*seed.DataFields.TS22) ;
                seed.DataFields.S22 = Es./(1-Nus.^2).*(seed.DataFields.TS22 + Nus.*seed.DataFields.TS11) ;
                seed.DataFields.S12 = Es./(1+Nus).*seed.DataFields.TS12 ;
                seed.DataFields.Seq = SIG ;
                
    %% Elasto-plastic with Isotropic hardening
        [S11,S22,S12] = Isotropic_Hardening(seed.DataFields.D11,seed.DataFields.D22,seed.DataFields.D12) ;
            seed.DataFields.CauchyStresses = 'Cauchy Stresses' ;
                seed.DataFields.S11 = S11 ;
                seed.DataFields.S22 = S22 ;
                seed.DataFields.S12 = S12 ;
                Smean = 1/3*(S11+S22) ;
                seed.DataFields.Seq = sqrt(3/2*((S11-Smean).^2 + (S22-Smean).^2 + Smean.^2 + 2*S12.^2)) ;
        
%% ADD ALL NEW FIELDS
    DATA = seed.DataFields ;
    DATA.PiolaKirchhoff2 = '2nd Piola-Kirchhoff Stresses' ;
        idet = (DATA.F11.*DATA.F22 - DATA.F21.*DATA.F12).^-1 ; 
        DATA.PK11 = idet.*(DATA.S11.*DATA.F22.^2 + DATA.S22.*DATA.F12.^2 - 2*DATA.S12.*DATA.F12.*DATA.F22) ;
        DATA.PK22 = idet.*(DATA.S11.*DATA.F21.^2 + DATA.S22.*DATA.F11.^2 - 2*DATA.S12.*DATA.F21.*DATA.F11) ;
        DATA.PK12 = idet.*(DATA.S12.*(DATA.F11.*DATA.F22 + DATA.F12.*DATA.F21) - DATA.S11.*DATA.F22.*DATA.F21 - DATA.S22.*DATA.F11.*DATA.F12) ;
    DATA.Power = 'Internal Power Density' ;
        DATA.pi_SxD = -(DATA.S11.*DATA.D11 + DATA.S22.*DATA.D22 + 2*DATA.S12.*DATA.D12) ;
        O = zeros(size(DATA.J,1),1) ;
        DATA.pi_PK2xLdot = -(DATA.PK11.*cat(3,O,diff(DATA.L11,1,3)) + DATA.PK22.*cat(3,O,diff(DATA.L22,1,3)) + 2.*DATA.PK12.*cat(3,O,diff(DATA.L12,1,3))) ;
    DATA.Energy = 'Strain Energy Density' ;
        DATA.w_SxD = cumsum(-DATA.pi_SxD,3) ;
        DATA.w_PK2xLdot = cumsum(-DATA.pi_PK2xLdot,3) ;
    seed.DataFields = DATA ;
    
%% INTERNAL POWER EQUALITY VERIFICATION

    clf ;
    seed = seed ;
    thickness = 0.8/1000 ;
    pixelsByMeters = 132000 ; % see below
    
    da = seed.integMat(seed.MovingPoints(:,:,1))*thickness/pixelsByMeters^2 ;
    
    Pi_SxD = -sum(DATA.pi_SxD.*DATA.J.*da(:),1,'omitnan') ;
    Pi_SeqxDeq = sum(DATA.Seq.*DATA.Deq.*DATA.J.*da(:),1,'omitnan') ;
    Pi_PK2xLdot = -sum(DATA.pi_PK2xLdot.*da(:),1,'omitnan') ;
    
    plot(Pi_SeqxDeq(:),'DisplayName','$\sigma_{eq} \, D_{eq}$') ;
    plot(Pi_SxD(:),'DisplayName','$\sigma_{ij} \, D_{ij}$') ;
    plot(Pi_PK2xLdot(:),'DisplayName','$\Pi_{ij} \, \dot L_{ij}$') ;
    lgd = legend ; 
    lgd.EdgeColor = 'k' ; 
    lgd.Orientation = 'horizontal' ;
    lgd.Location = 'best' ;
    
    % ADD THE POWER OF EXTERNAL FORCES
        % Boundary edges
            bndE = seed.Edges(seed.boundaryEdges,:) ;
        % Position of these edges's nodes
            x1e = reshape(seed.DataFields.x1(bndE,:,:),[],2,hd.nFrames) ;
            x2e = reshape(seed.DataFields.x2(bndE,:,:),[],2,hd.nFrames) ;
        % Normals (in the current config)
            normals = [-diff(x2e,1,2) diff(x1e,1,2)] ;
            dl = sqrt(sum(normals.^2,2)) ;
            normals = normals./dl ;
        % Vertical edges
            vertEdge = abs(normals(:,2,1)./normals(:,1,1))<.5 ;
        % Left & right edges
            x1Left = min(min(x1e(:,:,1),[],1)) ;
            x1Right = max(max(x1e(:,:,1),[],1)) ;
            L0 = x1Right-x1Left ;
            leftEdge = vertEdge & mean(x1e(:,:,1),2)<x1Left+0.1*L0 ;
            rightEdge = vertEdge & mean(x1e(:,:,1),2)>x1Right-0.1*L0 ;
        % Mean length
            L = mean(mean(x1e(rightEdge,:,:),2),1)-mean(mean(x1e(leftEdge,:,:),2),1) ;
            dL = [0 ; diff(L(:))] ;

        % Verify
%             clf
%             ii = 1 ;
%             axis ij ; axis equal ; axis off
%                 patch('Faces',seed.Elems,'Vertices',seed.MovingPoints(:,:,ii),'Facecolor','none','edgecolor','k')
%                 plot(x1e(vertEdge,:,ii)',x2e(vertEdge,:,ii)','k')
%                 plot(x1e(leftEdge,:,ii)',x2e(leftEdge,:,ii)','b')
%                 plot(x1e(rightEdge,:,ii)',x2e(rightEdge,:,ii)','r')
%             return ;

        F = hd.InputData ; F(isnan(F)) = 0 ;
        Pi_Fdl = F(:).*dL(:)/pixelsByMeters ;
        plot(Pi_Fdl,'DisplayName','$F \, \dot \ell$') ;
    

%% INITIALIZE THE FIGURES
    contourTag = 'J-integral Contours' ;
    valuesTag = 'J-integral Values' ;
    % Init the contour figure
        delete(findobj(0,'Name',contourTag)) ; 
        figContours = figure('Name',contourTag) ;
        refImg = hd.Images{1}(:,:,:,1) ;
        im = imagesc(repmat(refImg(:,:,1),[1 1 3])) ;
        axis tight
        axis equal
        axis ij
        axis off
        set(gca,'xtick',[],'ytick',[]) ; 
        [nI,nJ] = size(refImg) ;
    % Add a ruler line
        ruler = images.roi.Polyline ;
        ruler.Parent = gca ;
        ruler.Color = 'r' ;
        ruler.Position = [1/20 0.095 ; 1/20 0.865].*[nJ nI] ;
    % Init the values figure
        delete(findobj(0,'Name',valuesTag)) ; 
        figValues = figure('Name',valuesTag) ;
    % Position
        figValues.Position = [figValues.Position(1:2)+figValues.Position(3:4)/2.*[1 .5] figValues.Position(3:4)/2] ;
        figContours.Position = [figContours.Position(1:2)+figContours.Position(3:4)/2.*[0 .5] figContours.Position(3:4)/2] ;
        
%% SET THE INTEGRAL CONTOUR INTERACTIVELY AND COMPUTE THE INTEGRAL
    
    figure(figContours) ;
    crackPos = [0.5 0.80].*[nJ nI] ; % Crack tip position
    crackVec = [0 -1] ; % Crack direction
    zoneWidth = 0.30*nJ ;
    zoneHeight = 0.3*nI ;
    crackWidth = 0.05*nJ ;
    rulerLength = 0.007 ;
    
    nPts = 10000 ;
    
    pos = [crackPos(1)-crackWidth/2 crackPos(2)] ;
    pos(end+1,:) = [crackPos(1)-zoneWidth/2 crackPos(2)] ;
    pos(end+1,:) = [crackPos(1)-zoneWidth/2 crackPos(2)-zoneHeight] ;
    pos(end+1,:) = [crackPos(1)+zoneWidth/2 crackPos(2)-zoneHeight] ;
    pos(end+1,:) = [crackPos(1)+zoneWidth/2 crackPos(2)] ;
    pos(end+1,:) = [crackPos(1)+crackWidth/2 crackPos(2)] ;

    %delete(contour) ;
    contour = images.roi.Polyline ;
    contour.Parent = gca ;
    contour.Position = pos ;
    
    %contour.draw ; 


% COMPUTE THE J-INTEGRAL
    
    % Units
        pixelsByMeters = sqrt(sum(diff(ruler.Position,1,1).^2))/rulerLength ;
        
    % Integration points
        dL = sqrt(sum(diff(contour.Position,1,1).^2,2)) ;
        L0 = cumsum([0 ; dL]) ;
        s = linspace(0,L0(end),nPts) ; 
        Pts = interp1(L0,contour.Position,s) ;
        PtsMid = (Pts(1:end-1,:)+Pts(2:end,:))/2 ;
        dx = diff(Pts,1,1)/pixelsByMeters ;
        dl = sqrt(sum(dx.^2,2)) ;
        
    % Outgoing Normals
        normals = dx./dl*[0 1 ; -1 0]' ;
        switched = inpolygon(PtsMid(:,1)+normals(:,1),PtsMid(:,2)+normals(:,2),contour.Position(:,1),contour.Position(:,2)) ;
        normals(switched,:) = -normals(switched,:) ;
        
    % Interpolation matrix
        refPts = seed.MovingPoints(:,:,1) ;
        T = seed.interpMat(PtsMid,0,refPts) ;
        
    % Interpolated Data
        W = T*squeeze(cumsum(-seed.DataFields.pi_PK2xLdot,3)) ;
        du1_dx1 = T*squeeze(seed.DataFields.F11-1) ;
        du2_dx1 = T*squeeze(seed.DataFields.F21) ;
        du1_dx2 = T*squeeze(seed.DataFields.F12) ;
        du2_dx2 = T*squeeze(seed.DataFields.F22-1) ;
        PK11 = T*squeeze(seed.DataFields.PK11) ;
        PK22 = T*squeeze(seed.DataFields.PK22) ;
        PK12 = T*squeeze(seed.DataFields.PK12) ;
        t1 = normals(:,1).*PK11 + normals(:,2).*PK12 ;
        t2 = normals(:,1).*PK12 + normals(:,2).*PK22 ;
        
    % With respect to the crack direction
        du1_da = du1_dx1*crackVec(1) + du1_dx2*crackVec(2) ;
        du2_da = du2_dx1*crackVec(1) + du2_dx2*crackVec(2) ;
        na = sum(normals.*crackVec,2) ;
        
    % J integral
        j = W.*na - t1.*du1_da - t2.*du2_da ;
        J = sum(j.*dl) ;
        
    % Fracture toughness
        K = sqrt(J*180e9/(1-0.33^2)) ;
        max(K)
        
% PLOT THE RESULT
    figure(figValues) ;
    pl = plot(J/1000) ;
    %ii = 210 ; pl = plot(cumsum(dl)/sum(dl),j(:,ii)) ;
    contour.Color = pl.Color ;

    xlabel 'Image'
    ylabel 'J (kJ/m$^2$)'

    box on
    grid on
    
    
%% ADD UNLOADING FIELDS

    seed.DataFields.Deq = cat(3,...
                                                seed.DataFields.Leq(:,:,1)*0,...
                                                diff(seed.DataFields.Leq,1,3)) ;

    seed.DataFields.Dmax = seed.DataFields.Leq-max(seed.DataFields.Leq,[],3) ;
    
    seed.DataFields.logW = log10(seed.DataFields.W) ;
    
    
%% TRACE CONTOURS ON STRAIN ENERGY MAP
    
    if exist('figContours') && isvalid(figContours) ; close(figContours) ; end
    figContours = open('Jcontours.fig') ;
    if exist('figSE') && isvalid(figSE) ; close(figSE) ; end
    figSE = open('StrainEnergy.fig') ;
    
    axC = findobj(figContours,'type','axes') ;
    axSE = findobj(figSE,'type','axes') ;
    contours = findobj(axC.Children,'-property','InteractionsAllowed') ;
    
    iRef = 1 ;
    iCurrent = 282 ;
    refPts = seed.MovingPoints(:,:,iRef) ;
    pC = gobjects(0) ;
    for contour = contours(:)'
        dL = sqrt(sum(diff(contour.Position,1,1).^2,2)) ;
        L0 = cumsum([0 ; dL]) ;
        s = linspace(0,L0(end),nPts) ; 
        Pts = interp1(L0,contour.Position,s) ;
        interpMat = seed.interpMat(Pts,NaN,refPts) ;
        currPts = interpMat*seed.MovingPoints(:,:,iCurrent) ;
        pC(end+1) = plot(axSE,currPts(:,1),currPts(:,2),'linewidth',2,'color',contour.Color) ;
    end
    %set(pC,'color','k') ;
    
    
    
    
