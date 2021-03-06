%% J-INTEGRAL ESTIMATION ON DIC DATA


%% IMPORT THE SEED AND COMPUTE STRESSES AND STRAIN ENERGY
    seedNumber = 2 ;
    
    % Import hd struct
        global hd
        seed = hd.Seeds(seedNumber) ;
    % Re-compute the seed Data Fields on points
        seed.DataOnNodes = false ; true ;
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
        %E = 200e9 ; nu = .5 ; Sy = 427e6 ; K = 1175e6 ; n = .67 ; % MEAN ISOTROPIC
        E = 197e9 ; nu = .5 ; Sy = 419e6 ; K = 1136e6 ; n = .67 ; % ISN_0 (BUILD DIR.)
        %E = 201e9 ; nu = .5 ; Sy = 441e6 ; K = 1177e6 ; n = .68 ; % ISN_45 (OBLIQUE DIR.)
        %E = 201e9 ; nu = .5 ; Sy = 431e6 ; K = 1203e6 ; n = .63 ; % ISN_90 (PRINT DIR.)
        %E = 201e9 ; nu = .5 ; Sy = 493e6 ; K = 1065e6 ; n = .67 ; % ISN_90_LOADFIT (PRINT DIR.)
        [S11,S22,S33,S12,Dp11,Dp22,Dp33,Dp12,R] = Isotropic_Hardening(seed.DataFields.D11,seed.DataFields.D22,seed.DataFields.D12,E,nu,Sy,K,n) ;
            seed.DataFields.ElasticStrains = 'Elastic Strains' ;
                seed.DataFields.De11 = seed.DataFields.D11-Dp11 ;
                seed.DataFields.De22 = seed.DataFields.D22-Dp22 ;
                seed.DataFields.De33 = seed.DataFields.D33-Dp33 ;
                seed.DataFields.De12 = seed.DataFields.D12-Dp12 ;
                seed.DataFields.TSe11 = cumsum(seed.DataFields.De11,3) ;
                seed.DataFields.TSe22 = cumsum(seed.DataFields.De22,3) ;
                seed.DataFields.TSe33 = cumsum(seed.DataFields.De33,3) ;
                seed.DataFields.TSe12 = cumsum(seed.DataFields.De12,3) ;
            seed.DataFields.PlasticStrains = 'Plastic Strains' ;
                seed.DataFields.Dp11 = Dp11 ;
                seed.DataFields.Dp22 = Dp22 ;
                seed.DataFields.Dp33 = Dp33 ;
                seed.DataFields.Dp12 = Dp12 ;
                Dmean = 1/3*(seed.DataFields.Dp11+seed.DataFields.Dp22+seed.DataFields.Dp33) ;
                seed.DataFields.dP = sqrt(2/3*((seed.DataFields.Dp11-Dmean).^2 + (seed.DataFields.Dp22-Dmean).^2 + (seed.DataFields.Dp33-Dmean).^2 + 2*(seed.DataFields.Dp12-Dmean).^2)) ;
                seed.DataFields.P = cumsum(seed.DataFields.dP,3) ;
                seed.DataFields.TSp11 = cumsum(Dp11,3) ;
                seed.DataFields.TSp22 = cumsum(Dp22,3) ;
                seed.DataFields.TSp33 = cumsum(Dp33,3) ;
                seed.DataFields.TSp12 = cumsum(Dp12,3) ;
            seed.DataFields.CauchyStresses = 'Cauchy Stresses' ;
                seed.DataFields.S11 = S11 ;
                seed.DataFields.S22 = S22 ;
                seed.DataFields.S12 = S12 ;
                Smean = 1/3*(S11+S22) ;
                seed.DataFields.Seq = sqrt(3/2*((S11-Smean).^2 + (S22-Smean).^2 + Smean.^2 + 2*S12.^2)) ;
                seed.DataFields.R = R ;
                seed.DataFields.f = (seed.DataFields.Seq-seed.DataFields.R)./seed.DataFields.R ;
                seed.DataFields.Dplast = seed.DataFields.R.*seed.DataFields.dP ;
                
        %% STRESS REGULARISATION
            mu = E/2/(1+nu) ;
            nPts = size(seed.MovingPoints,1) ;
            nElems = size(seed.Elems,1) ;
            lambda = NaN(nPts,2,hd.nFrames) ;
            regularize = 'total' ; 'rate' ; % regularize the total or rate of elastic strains
            % Boundary points (on which lambda is zero)
                bndPts = reshape(seed.Edges(seed.boundaryEdges,:),[],1) ;
                bndPts = ismember(1:nPts,bndPts) ;
            % Elastic strain rate
                De = [  seed.DataFields.D11(:,:)-seed.DataFields.Dp11(:,:) ; ...
                        seed.DataFields.D12(:,:)-seed.DataFields.Dp12(:,:) ; ...
                        seed.DataFields.D12(:,:)-seed.DataFields.Dp12(:,:) ; ...
                        seed.DataFields.D22(:,:)-seed.DataFields.Dp22(:,:) ; ...
                        seed.DataFields.D33(:,:)-seed.DataFields.Dp33(:,:) ; ...
                        ] ;
            % Corrector
                Dc = De*NaN ;
            % REGULARIZATION
                wtbr = waitbar(0,'Computing Regularization...') ;
                for ii = 1:hd.nFrames
                    % Data to regularize
                        switch regularize
                            case 'rate'
                                ee = De(:,ii)-kron([1;0;0;1;1],De(end-nElems+1:end,ii)) ;
                            case 'total'
                                ee = sum(De(:,1:ii)-kron([1;0;0;1;1],De(end-nElems+1:end,1:ii)),2) ;
                        end
                    % Gradient matrices
                        [D1,D2] = seed.diffMat(seed.MovingPoints(:,:,ii)) ;
                        D1(:,bndPts) = [] ; D2(:,bndPts) = [] ;
                        O = D1*0 ;
                    % Differential operators
                        B = [D1 O ; O D1 ; D2 O ; O D2 ; O O] ;
                        Bt = [D1 O ; D2 O ; O D1 ; O D2 ; O O] ;
                        Bdiv = [O O ; O O ; O O ; O O ; D1 D2] ;
                    % Problem assembly
                        K = B'*(B+3*Bt) ;
                        f = -2*B'*ee ;
                    % Solve for lambda
                        ll = K\f ;
                        lambda(~bndPts,:,ii) = reshape(ll,[],2) ;
                        lambda(bndPts,:,ii) = 0 ;
                    % Compute the elastic correction
                        Dc(:,ii) = (1/2*(B + Bt) - Bdiv)*ll ;
                        if strcmp(regularize,'total') && ii>1
                            %Dc(:,ii) = Dc(:,ii) - sum(Dc(:,1:ii-1),2) ;
                        end
                        De(:,ii) = De(:,ii) + Dc(:,ii) ;
                    % Waitbar
                        wtbr = waitbar(ii/hd.nFrames,wtbr,['Computing Regularization...(' num2str(ii) '/' num2str(hd.nFrames) ')']) ;
                end
                delete(wtbr) ;
            % NEW FIELDS
                seed.DataFields.Regularized = 'Stress Regularization' ;
                    seed.DataFields.LAMBDA1 = reshape(lambda(:,1,:),nPts,1,hd.nFrames) ;
                    seed.DataFields.LAMBDA2 = reshape(lambda(:,2,:),nPts,1,hd.nFrames) ;
                    seed.DataFields.LAMBDA = sqrt(seed.DataFields.LAMBDA1.^2 + seed.DataFields.LAMBDA2.^2) ;
                    seed.DataFields.Dc11 = reshape(Dc(1:nElems,:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Dc12 = reshape(Dc(nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Dc22 = reshape(Dc(3*nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Dc33 = reshape(Dc(4*nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Dec11 = reshape(De(1:nElems,:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Dec12 = reshape(De(nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Dec22 = reshape(De(3*nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Dec33 = reshape(De(4*nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    seed.DataFields.TSec11 = cumsum(seed.DataFields.Dec11,3) ;
                    seed.DataFields.TSec12 = cumsum(seed.DataFields.Dec12,3) ;
                    seed.DataFields.TSec22 = cumsum(seed.DataFields.Dec22,3) ;
                    seed.DataFields.TSec33 = cumsum(seed.DataFields.Dec33,3) ;
                    seed.DataFields.dTSe11 = seed.DataFields.TSec11-seed.DataFields.TSe11 ;
                    seed.DataFields.dTSe12 = seed.DataFields.TSec12-seed.DataFields.TSe12 ;
                    seed.DataFields.dTSe22 = seed.DataFields.TSec22-seed.DataFields.TSe22 ;
                    seed.DataFields.dTSe33 = seed.DataFields.TSec33-seed.DataFields.TSe33 ;
                    TSe = cumsum(De,2) ;
                    Sr = 2*mu*(TSe-kron([1;0;0;1;1],TSe(end-nElems+1:end,:))) ;
                    seed.DataFields.Sr11 = reshape(Sr(0*nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Sr22 = reshape(Sr(3*nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    seed.DataFields.Sr12 = reshape(Sr(1*nElems+(1:nElems),:),nElems,1,hd.nFrames) ;
                    SrMean = 1/3*(seed.DataFields.Sr11+seed.DataFields.Sr22) ;
                    seed.DataFields.Sreq = sqrt(3/2*((seed.DataFields.Sr11-SrMean).^2 + (seed.DataFields.Sr22-SrMean).^2 + SrMean.^2 + 2*seed.DataFields.Sr12.^2)) ;
                    seed.DataFields.fr = (seed.DataFields.Sreq - seed.DataFields.R)./seed.DataFields.R ;
                    seed.DataFields.dS11 = seed.DataFields.Sr11-seed.DataFields.S11 ;
                    seed.DataFields.dS22 = seed.DataFields.Sr22-seed.DataFields.S22 ;
                    seed.DataFields.dS12 = seed.DataFields.Sr12-seed.DataFields.S12 ;
                    seed.DataFields.Drplast = seed.DataFields.Sreq.*seed.DataFields.dP ;
                
                
%% CHECK DIVERGENCE OF STRESSES
    [D1,D2] = seed.diffMat(seed.MovingPoints(:,:,1)) ;
    seed.DataFields.Divergence = 'Stress Divergence' ;
        pixelsByMeters = 135000 ; % see below
        seed.DataFields.f1 = -reshape(D1*S11(:,:) + D2*S12(:,:),[],1,hd.nFrames)*pixelsByMeters ;
        seed.DataFields.f2 = -reshape(D1*S12(:,:) + D2*S22(:,:),[],1,hd.nFrames)*pixelsByMeters ;
        
        
%% CHECK RESULTANT LOAD

    thickness = .6/1000 ;
    pixelsByMeters = 135000 ; % see below
    S11 = seed.DataFields.Sr11 ; 
    S22 = seed.DataFields.Sr22 ; 
    S12 = seed.DataFields.Sr12 ;
    
    % COMPUTE THE LOAD
        % Boundary edges
            bndE = seed.Edges(seed.boundaryEdges,:) ;
        % Position of these edges's nodes
            x1e = reshape(seed.DataFields.x1(bndE,:,:),[],2,hd.nFrames) ;
            x2e = reshape(seed.DataFields.x2(bndE,:,:),[],2,hd.nFrames) ;
        % Normals (in the current config)
            normals = [-diff(x2e,1,2) diff(x1e,1,2)] ;
            dl = sqrt(sum(normals.^2,2)) ;
            normals = normals./dl ;
        % Outgoing normals
            ps = [mean(x1e(:,:,1),2) mean(x2e(:,:,1),2)] + 0.1*normals(:,:,1) ;
            outside = any(isnan(seed.interpMat(ps,NaN,seed.MovingPoints(:,:,1))),2) ;
            normals(~outside,:,:) = -normals(~outside,:,:) ;
        % Vertical edges
            vertEdge = abs(normals(:,2,1)./normals(:,1,1))<.5 ;
        % Left & right edges
            x1Left = min(min(x1e(:,:,1),[],1)) ;
            x1Right = max(max(x1e(:,:,1),[],1)) ;
            L0 = x1Right-x1Left ;
            leftEdge = vertEdge & mean(x1e(:,:,1),2)<x1Left+0.1*L0 ;
            rightEdge = vertEdge & mean(x1e(:,:,1),2)>x1Right-0.1*L0 ;
        % Corresponding elements
            [~,elem2edg] = seed.getEdges ;
            elem2bndEdg = elem2edg(seed.boundaryEdges,:) ;
            elems2leftEdg = elem2bndEdg(leftEdge,:) ;
            elems2rightEdg = elem2bndEdg(rightEdge,:) ;
            leftElems = any(elems2leftEdg,1) ;
            rightElems = any(elems2rightEdg,1) ;
        % Resultant
            n1 = normals(:,1,:) ; n2 = normals(:,2,:) ;
            da = dl/pixelsByMeters*thickness ;
            Tl1 = da(leftEdge).*[n1(leftEdge,:).*(elems2leftEdg*S11(:,:)) + n2(leftEdge,:).*(elems2leftEdg*S12(:,:))] ;
            Tl2 = da(leftEdge).*[n1(leftEdge,:).*(elems2leftEdg*S12(:,:)) + n2(leftEdge,:).*(elems2leftEdg*S22(:,:))] ;
            Tr1 = da(rightEdge).*[n1(rightEdge,:).*(elems2rightEdg*S11(:,:)) + n2(rightEdge,:).*(elems2rightEdg*S12(:,:))] ;
            Tr2 = da(rightEdge).*[n1(rightEdge,:).*(elems2rightEdg*S12(:,:)) + n2(rightEdge,:).*(elems2rightEdg*S22(:,:))] ;
            Fl1 = sum(Tl1,1) ;
            Fl2 = sum(Tl2,1) ;
            Fr1 = sum(Tr1,1) ;
            Fr2 = sum(Tr2,1) ;

        % Verify
%             clf
%             ii = 260 ;
%             axis ij ; axis equal ; axis off
%                 patch('Faces',seed.Elems,'Vertices',seed.MovingPoints(:,:,ii),'Facecolor','none','edgecolor','k')
%                 plot(x1e(vertEdge,:,ii)',x2e(vertEdge,:,ii)','k')
%                 plot(x1e(leftEdge,:,ii)',x2e(leftEdge,:,ii)','b')
%                 plot(x1e(rightEdge,:,ii)',x2e(rightEdge,:,ii)','r')
%                 patch('Vertices',seed.MovingPoints(:,:,ii),'Faces',seed.Elems(leftElems,:),'facecolor','b') ;
%                 patch('Vertices',seed.MovingPoints(:,:,ii),'Faces',seed.Elems(rightElems,:),'facecolor','r') ;
% %                 quiver(mean(x1e(leftEdge,:,ii),2),mean(x2e(leftEdge,:,ii),2),normals(leftEdge,1,ii),normals(leftEdge,2,ii),'b') ;
% %                 quiver(mean(x1e(rightEdge,:,ii),2),mean(x2e(rightEdge,:,ii),2),normals(rightEdge,1,ii),normals(rightEdge,2,ii),'r') ;
%                 quiver(mean(x1e(leftEdge,:,ii),2),mean(x2e(leftEdge,:,ii),2),Tl1(:,ii),Tl2(:,ii),'b') ;
%                 quiver(mean(x1e(rightEdge,:,ii),2),mean(x2e(rightEdge,:,ii),2),Tr1(:,ii),Tr2(:,ii),'r') ;
%             return ;

    
%     clf
%         plot(hd.InputData,'k','DisplayName','measured') ;
%         plot(-Fl1,'DisplayName','$-F_1$ left') ;
%         plot(Fr1,'DisplayName','$F_1$ right') ;
%         plot(Fl2,'DisplayName','$F_2$ left') ;
%         plot(Fr2,'DisplayName','$F_2$ right') ;
%         legend('edgecolor','k')
    
        %plot(-Fl1,'DisplayName','$-F_1$ left (iden+reg)') ;%+reg
        %plot(Fr1,'DisplayName','$F_1$ right (iden+reg)') ;
        plot((Fr1-Fl1)/2,'DisplayName','$F$ (optim)') ;
        set(findobj(gca,'type','line'),'linewidth',1) ;
    
        
%% ADD ALL NEW FIELDS
    DATA = seed.DataFields ;
    S11 = DATA.Sr11 ; S22 = DATA.Sr22 ; S12 = DATA.Sr12 ;
    DATA.PiolaKirchhoff2 = '2nd Piola-Kirchhoff Stresses' ;
        idet = (DATA.F11.*DATA.F22 - DATA.F21.*DATA.F12).^-1 ; 
        DATA.PKII11 = idet.*(S11.*DATA.F22.^2 + S22.*DATA.F12.^2 - 2*S12.*DATA.F12.*DATA.F22) ;
        DATA.PKII22 = idet.*(S11.*DATA.F21.^2 + S22.*DATA.F11.^2 - 2*S12.*DATA.F21.*DATA.F11) ;
        DATA.PKII12 = idet.*(S12.*(DATA.F11.*DATA.F22 + DATA.F12.*DATA.F21) - S11.*DATA.F22.*DATA.F21 - S22.*DATA.F11.*DATA.F12) ;
    DATA.PiolaKirchhoff1 = '1st Piola-Kirchhoff Stresses' ;
        DATA.PKI11 = DATA.F11.*DATA.PKII11 + DATA.F12.*DATA.PKII12 ;
        DATA.PKI12 = DATA.F11.*DATA.PKII12 + DATA.F12.*DATA.PKII22 ;
        DATA.PKI21 = DATA.F21.*DATA.PKII11 + DATA.F22.*DATA.PKII12 ;
        DATA.PKI22 = DATA.F21.*DATA.PKII12 + DATA.F22.*DATA.PKII22 ;
    DATA.Power = 'Internal Power Density' ;
        DATA.pi_SxD = -(S11.*DATA.D11 + S22.*DATA.D22 + 2*S12.*DATA.D12) ;
        O = zeros(size(DATA.J,1),1) ;
        DATA.pi_PK1xFdot = -(DATA.PKI11.*cat(3,O,diff(DATA.F11,1,3)) + DATA.PKI22.*cat(3,O,diff(DATA.F22,1,3)) + DATA.PKI12.*cat(3,O,diff(DATA.F12,1,3)) + DATA.PKI21.*cat(3,O,diff(DATA.F21,1,3))) ;
        DATA.pi_PK2xLdot = -(DATA.PKII11.*cat(3,O,diff(DATA.L11,1,3)) + DATA.PKII22.*cat(3,O,diff(DATA.L22,1,3)) + 2.*DATA.PKII12.*cat(3,O,diff(DATA.L12,1,3))) ;
    DATA.Energy = 'Strain Energy Density' ;
        DATA.w_SxD = cumsum(-DATA.pi_SxD,3) ;
        DATA.w_PK2xLdot = cumsum(-DATA.pi_PK2xLdot,3) ;
    seed.DataFields = DATA ;
    
    
%% INTERNAL POWER EQUALITY VERIFICATION

    clf ;
    seed = seed ;
    thickness = 0.6/1000 ;
    pixelsByMeters = 135000 ; % see below
    DATA = seed.DataFields ;
    
    da = seed.integMat(seed.MovingPoints(:,:,1))*thickness/pixelsByMeters^2 ;
    
    Pi_SxD = -sum(DATA.pi_SxD.*DATA.J.*da(:),1,'omitnan') ;
    Pi_SeqxDeq = sum(DATA.Sreq.*DATA.Deq.*DATA.J.*da(:),1,'omitnan') ;
    Pi_PK1xFdot = -sum(DATA.pi_PK1xFdot.*da(:),1,'omitnan') ;
    Pi_PK2xLdot = -sum(DATA.pi_PK2xLdot.*da(:),1,'omitnan') ;
    Pi_Plast = sum(DATA.Drplast.*DATA.J.*da(:),1,'omitnan') ;
    
    %plot(Pi_SeqxDeq(:),'DisplayName','$\mathcal{P}_i = -\int \sigma_{eq} \, D_{eq} \, \mathrm{d}\Omega$') ;
    plot(Pi_SxD(:),'DisplayName','$\mathcal{P}_i = -\int \sigma_{ij} \, D_{ij} \, \mathrm{d}\Omega$') ;
    plot(Pi_PK1xFdot(:),'DisplayName','$\mathcal{P}_i = -\int \Pi^1_{ij} \, \dot F_{ij} \, \mathrm{d}\Omega_0$') ;
    plot(Pi_PK2xLdot(:),'DisplayName','$\mathcal{P}_i = -\int \Pi^2_{ij} \, \dot L_{ij} \, \mathrm{d}\Omega_0$') ;
    plot(Pi_Plast(:),'DisplayName','$\mathcal{D}_\mathrm{plast} = -\int \sigma_{eq} \, \dot p \, \mathrm{d}\Omega$') ;
    lgd = legend ; 
    lgd.EdgeColor = 'k' ; 
    lgd.Orientation = 'horizontal' ;
    lgd.Location = 'best' ;
    lgd.Orientation = 'vertical' ;
    
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
        % Power
            F = hd.InputData ; F(isnan(F)) = 0 ;
            Pi_Fdl = F(:).*dL(:)/pixelsByMeters ;
            plot(Pi_Fdl,'DisplayName','$\mathcal{P}_e = F \, \dot \ell$') ;

        % Verify
%             clf
%             ii = 1 ;
%             axis ij ; axis equal ; axis off
%                 patch('Faces',seed.Elems,'Vertices',seed.MovingPoints(:,:,ii),'Facecolor','none','edgecolor','k')
%                 plot(x1e(vertEdge,:,ii)',x2e(vertEdge,:,ii)','k')
%                 plot(x1e(leftEdge,:,ii)',x2e(leftEdge,:,ii)','b')
%                 plot(x1e(rightEdge,:,ii)',x2e(rightEdge,:,ii)','r')
    

%% INITIALIZE THE FIGURES
    contourTag = 'J-integral Contours' ;
    valuesTag = 'J-integral Values' ;
    % Init the contour figure
        figContours = findobj(0,'Name',contourTag) ;
        if isempty(figContours)
            figContours = figure('Name',contourTag) ; 
            figContours.Position = [figContours.Position(1:2)+figContours.Position(3:4)/2.*[0 .5] figContours.Position(3:4)/2] ;
        end
        clf(figContours) ;
        refImg = hd.Images{1}{1} ;
        im = imagesc(repmat(refImg(:,:,1),[1 1 3])) ;
        mesh = patch('Vertices',seed.MovingPoints(:,:,1),'Faces',seed.Elems,'Facecolor','w','edgecolor','none','Facealpha',0.5) ;
        axis tight
        axis equal
        axis ij
        axis off
        set(gca,'xtick',[],'ytick',[]) ; 
        [nI,nJ] = size(refImg) ;
    % Add a ruler line
        ruler = images.roi.Line ;
        ruler.Parent = gca ;
        ruler.Color = 'r' ;
        ruler.Position = [1415 1130 ; 1550 1130] ; %[1/20 0.095 ; 1/20 0.865].*[nJ nI] ;
    % Init the values figure
        figValues = findobj(0,'Name',valuesTag) ;
        if isempty(figValues)
            figValues = figure('Name',valuesTag) ; 
            figValues.Position = [figValues.Position(1:2)+figValues.Position(3:4)/2.*[1 .5] figValues.Position(3:4)/2] ;
        end
        clf(figValues) ;
        
%% SET THE INTEGRAL CONTOUR INTERACTIVELY AND COMPUTE THE INTEGRAL
    
    figure(figContours) ;
    crackPos = [0.5 0.825].*[nJ nI] ; % Crack tip position [0.84 0.835 0.83 0.825]
    crackVec = [0 -1] ; % Crack direction
    zoneWidth = 0.25*nJ ; % [0.7 0.55 0.4 0.25]
    zoneHeight = 0.69*nI ; % [0.72 0.71 0.70 0.69]
    crackWidth = 0.05*nJ ;
    rulerLength = 0.001 ;
    iRef = 1 ; % Reference image ;
    configuration = 'current' ; % 'reference' is not working
    
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
        
    % Integration points (in the reference config.)
        dL = sqrt(sum(diff(pos,1,1).^2,2)) ;
        L0 = cumsum([0 ; dL]) ;
        s = linspace(0,L0(end),nPts) ; 
        refPts = interp1(L0,pos,s) ;
        
    % Integration points (current config)
        Tn = seed.interpMat(refPts,NaN,seed.MovingPoints(:,:,iRef)) ;
        switch configuration
            case 'current'
                Pts = reshape(Tn*reshape(seed.MovingPoints,[],2*hd.nFrames),nPts,2,hd.nFrames) ;
            case 'reference'
                Pts = Tn*seed.MovingPoints(:,:,iRef) ;
        end
        midPts = (Pts(1:end-1,:,:)+Pts(2:end,:,:))/2 ;
        
    % Outgoing Normals
        dx = diff(Pts,1,1)/pixelsByMeters ;
        dl = sqrt(sum(dx.^2,2)) ;
        normals = [dx(:,2,:) , -dx(:,1,:)]./dl ;
        switched = inpolygon(midPts(:,1,iRef)+normals(:,1,iRef),midPts(:,2,iRef)+normals(:,2,iRef),Pts(:,1,iRef),Pts(:,2,iRef)) ;
        normals(switched,:,:) = -normals(switched,:,:) ;
        n1 = squeeze(normals(:,1,:)) ; n2 = squeeze(normals(:,2,:)) ;
        
    % Verify
%         clf ;
%             mesh = patch('Vertices',seed.Points,'Faces',seed.Elems,'facecolor','none','edgecolor','k') ;
%             line = plot(Pts(:,1,1),Pts(:,2,1),'k') ;
%             qui = quiver(midPts(:,1,1),midPts(:,2,1),n1(:,1),n2(:,1)) ;
%             axis equal, axis ij, axis tight, axis off
%             for ii = 1:hd.nFrames
%                 mesh.Vertices = seed.MovingPoints(:,:,ii) ;
%                 line.XData = Pts(:,1,ii) ;
%                 line.YData = Pts(:,2,ii) ;
%                 qui.XData = midPts(:,1,ii) ; qui.YData = midPts(:,2,ii) ;
%                 qui.UData = n1(:,ii) ; qui.VData = n2(:,ii) ;
%                 drawnow ;
%             end
%         return
        
    % Interpolated Data
        % Interpolation matrix (mean over edges)
            Te = seed.interpMat((refPts(1:end-1,:)+refPts(2:end,:))/2,NaN,seed.MovingPoints(:,:,iRef)) ;
            if size(seed.DataFields.J,1)==size(seed.Elems,1) ; Te = Te*seed.elem2nod ; end
        % Inverse transformation
            iF11 = seed.DataFields.F22./seed.DataFields.J ;
            iF22 = seed.DataFields.F11./seed.DataFields.J ;
            iF21 = -seed.DataFields.F21./seed.DataFields.J ;
            iF12 = -seed.DataFields.F12./seed.DataFields.J ;
        switch configuration
            case 'current'
                W = cumsum(-seed.DataFields.pi_SxD,3) ;
%                 du1_dx1 = 1-iF11 ;
%                 du2_dx1 = -iF21 ;
%                 du1_dx2 = -iF12 ;
%                 du2_dx2 = 1-iF22 ;
                du1_dx1 = seed.DataFields.F11-1 ;
                du2_dx1 = seed.DataFields.F21 ;
                du1_dx2 = seed.DataFields.F12 ;
                du2_dx2 = seed.DataFields.F22-1 ;
                S11 = seed.DataFields.Sr11 ;
                S22 = seed.DataFields.Sr22 ;
                S12 = seed.DataFields.Sr12 ;
                S21 = seed.DataFields.Sr12 ;
            case 'reference'
                W = cumsum(-seed.DataFields.pi_PK2xLdot,3) ;
%                 du1_dx1 = seed.DataFields.F11.*iF11 + seed.DataFields.F21.*iF21 ;
%                 du2_dx1 = seed.DataFields.F12.*iF11 + seed.DataFields.F22.*iF21 ;
%                 du1_dx2 = seed.DataFields.F11.*iF12 + seed.DataFields.F21.*iF22 ;
%                 du2_dx2 = seed.DataFields.F12.*iF12 + seed.DataFields.F22.*iF22 ;
                du1_dx1 = seed.DataFields.F11-1 ;
                du2_dx1 = seed.DataFields.F21 ;
                du1_dx2 = seed.DataFields.F12 ;
                du2_dx2 = seed.DataFields.F22-1 ;
%                 S11 = -seed.DataFields.PKII11 ;
%                 S22 = -seed.DataFields.PKII22 ;
%                 S12 = -seed.DataFields.PKII12 ;
%                 S21 = -seed.DataFields.PKII12 ;
                S11 = seed.DataFields.PKI11 ;
                S22 = seed.DataFields.PKI22 ;
                S12 = seed.DataFields.PKI21 ;
                S21 = seed.DataFields.PKI12 ;
        end
        
    % Stress vector
        t1 = n1.*(Te*squeeze(S11)) + n2.*(Te*squeeze(S21)) ;
        t2 = n1.*(Te*squeeze(S12)) + n2.*(Te*squeeze(S22)) ;
        
    % Projection in the crack direction
        du1_da = (Te*squeeze(du1_dx1))*crackVec(1) + (Te*squeeze(du1_dx2))*crackVec(2) ;
        du2_da = (Te*squeeze(du2_dx1))*crackVec(1) + (Te*squeeze(du2_dx2))*crackVec(2) ;
        na = squeeze(sum(normals.*crackVec,2)) ;
        
    % J integral
        j = (Te*squeeze(W)).*na - t1.*du1_da - t2.*du2_da ;
        J = sum(j.*dl(:,:),1,'omitnan') ;
        J(all(isnan(j),1)) = NaN ;
        
    % Fracture toughness
        K = sqrt(J*E/(1-nu^2)) ;
        max(K)
        K(210)
        
% PLOT THE RESULT
    figure(figValues) ;
    pl = plot(J/1000) ;
    %ii = 210 ; pl = plot(cumsum(dl)/sum(dl),j(:,ii)) ;
    contour.Color = pl.Color ;

    xlabel 'Image'
    ylabel 'J (kJ/m$^2$)'

    box on
    grid on
    
    
%% FOLLOW THE CRACK TIP
    iRef = 1 ;
    u1 = seed.DataFields.u1(:,:) ;
    
    Rmax = 100 ; % Query zone radius
    xc = [800 450] ; % Crack tip position
    t0 = pi/2 ; % Crack direction (angle)
    dt = 90*pi/180 ; % Crack opening angle
    uh = 0 ; % Homogeneous displacement
    kappa = (3-nu)/(1+nu) ;
    
    p0 = [xc dt uh] ; % parameters
    xc = @(p) p(1:2) ;
    dt = @(p) p(3) ;
    uh = @(p) p(4) ;
    
    x = seed.MovingPoints(:,:,iRef) ;
    
    xp = @(p) (x-xc(p))*[cos(t0) sin(t0) ; -sin(t0) cos(t0)] ;
    r = @(p) sqrt(sum(xp(p).^2,2)) ;
    t = @(p) angle(xp(p)*[1 ; 1i]) ;
    valid = @(p) r(p)<Rmax ;
    dist = @(p,q) norm((u(p)-q).*valid(p)) ;
    
    % u = u0 + t/2/pi.*sin(dt).*r ;
    u = @(p)  uh(p) + sin(dt(p))*sqrt(r(p)).*((2*kappa+1).*sin(t(p)/2)-sin(3*t(p)/2)) ;
    
    u0 = u(p0) ;

    cla ; 
    ax = gca ;
        axis equal
        axis ij
    mesh = patch('Faces',seed.Elems,'facecolor','none','edgecolor','k') ;
    pts = plot3(x(valid(p0),1),x(valid(p0),2),u0(valid(p0)),'.r') ;
    for ii = 240%1:hd.nFrames
        q = u1(:,ii) ;
        p = fminsearch(@(p)dist(p,q),p0) ;
        mesh.Vertices = [x q] ;
        mesh.Vertices(~valid(p),:) = NaN ;
        pts.XData = x(valid(p),1) ;
        pts.YData = x(valid(p),2) ;
        ui = u(p) ; pts.ZData = ui(valid(p)) ;
        drawnow ;
    end
        
    
    
%% TRACE CONTOURS ON STRAIN ENERGY MAP
    
    if exist('figContours') && isvalid(figContours) ; close(figContours) ; end
    figContours = open('Jcontours.fig') ;
    if exist('figSE') && isvalid(figSE) ; close(figSE) ; end
    figSE = open('RegOnsetStrainEnergy.fig') ;
    cornerID = 4 ; % For labels
    iRef = 1 ;
    iCurrent = 209 ; 282 ;
    
    axC = findobj(figContours,'type','axes') ;
    axSE = findobj(figSE,'type','axes') ;
    contours = findobj(axC.Children,'-property','InteractionsAllowed') ;
    
    refPts = seed.MovingPoints(:,:,iRef) ;
    pC = gobjects(0) ;
    txt = gobjects(0) ;
    for contour = contours(:)'
        if size(contour.Position,1)<3 ; continue ; end
        dL = sqrt(sum(diff(contour.Position,1,1).^2,2)) ;
        L0 = cumsum([0 ; dL]) ;
        s = linspace(0,L0(end),nPts) ; 
        refPts = interp1(L0,contour.Position,s) ;
        interpMat = seed.interpMat(refPts,NaN,refPts) ;
        currPts = interpMat*seed.MovingPoints(:,:,iCurrent) ;
        pC(end+1) = plot(axSE,currPts(:,1),currPts(:,2),'linewidth',1,'color',contour.Color) ;
        corner = seed.interpMat(contour.Position(cornerID,:),0,refPts)*seed.MovingPoints(:,:,iCurrent) ;
        txt(end+1) = text(axSE,corner(1),corner(2),['' char(length(txt)+97) '.'],'color',contour.Color) ;
    end
    set(pC,'color','k') ;
    set(txt,'color','k','verticalalignment','top','horizontalalignment','right','fontname','consolas','interpreter','tex') ;
    
    close(figContours) ;
    
    
    
    
