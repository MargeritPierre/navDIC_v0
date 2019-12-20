%% ANALYSIS OF A PERIODIC STRUCTURE

%% CREATE THE NECESSARY MASKS

    global hd
    clearvars -except hd

    macroSeed = hd.Seeds(1) ;
    microSeed = hd.Seeds(2) ;
    macroPositionFrom = microSeed ; macroSeed ;
    nodesByQuad = 4 ; % 4 or 8 
    nPx = 40 ; nPy = nPx ;
    refFrame = 1 ;
    IMG = @(ii)repmat(hd.Images{1}(:,:,:,ii),[1 1 3]) ;

    indQUADS = logical([0 1 1 1 0]'*[0 0 1 1 1 1 0 0]) ;
    nQuads = sum(indQUADS(:)) ;
    QUADS = macroSeed.Quadrangles(indQUADS(:),:) ;
    
    % Create the original grid of element coordinates
        [Ex,Ey] = meshgrid(linspace(0,1,nPx),linspace(0,1,nPy)) ;
        ex = Ex(:)' ; ey = Ey(:)' ;
        
    % Position of the macro quad nodes
        switch macroPositionFrom
            case microSeed
                P0 = macroSeed.MovingPoints(QUADS,:,refFrame) ;
                Pq = microSeed.interpMat(P0)*reshape(microSeed.MovingPoints,size(microSeed.MovingPoints,1),[]) ;
                Pq = reshape(Pq,[],4,2,hd.nFrames) ;
            case macroSeed
                Pq = reshape(macroSeed.MovingPoints(QUADS,:,:),[],4,2,hd.nFrames) ;
        end
        
    % Quad mean dimensions (in the reference config)
        Lq = squeeze(mean(range(Pq(:,:,:,refFrame),2))) ;
        
    % Macro Position interpolated in the image coordinates
        PPPi = Pq(:,1,:,:).*(1-ex).*(1-ey) ...
                + Pq(:,2,:,:).*(ex).*(1-ey) ...
                + Pq(:,3,:,:).*(ex).*(ey) ...
                + Pq(:,4,:,:).*(1-ex).*(ey) ;
    
    % Create the interpolation matrix
        PPPi0 = reshape(PPPi(:,:,:,refFrame),[],2) ; % Macro position in the ref. config.
        interpMat = microSeed.interpMat(PPPi0,NaN) ;
        
    % Points INSIDE the micro seed mesh
        vMASK = ~any(isnan(interpMat),2) ; % Pin = Pe(MASK,:) ;
        
    % Take the same points for all quads
        qMASK = reshape(vMASK,nQuads,nPx*nPy) ; % Individual quad mask
        tMASK = all(qMASK,1) ; % Points valid in all quads
        cMASK = qMASK & tMASK ; % Common masks
        nCommonPts = sum(tMASK(:)) ;
        
    % Display
        clf ;
            im = imagesc(IMG(1)) ;
            axis ij
            axis tight
            axis equal
            plIn = plot(PPPi0(vMASK,1),PPPi0(vMASK,2),'.b','markersize',5) ;
            plCommonIn = plot(PPPi0(cMASK,1),PPPi0(cMASK,2),'.c','markersize',5) ;
            plQ = plot(Pq(:,:,1),Pq(:,:,2),'.r','markersize',20) ;
            %plOut = plot(PPPi0(~vMASK,1),PPPi0(~vMASK,2),'.g','markersize',1) ;
            
            
%% INTERPOLATION OF MICRO DISPLACEMENTS

    % Retrieve MICRO Positions of the interpolated points
        interpMat(~cMASK(:),1) = NaN ;
        PPP = reshape(interpMat*reshape(microSeed.MovingPoints,size(interpMat,2),[]),[],2,hd.nFrames) ;
        UUU = PPP-PPP(:,:,refFrame) ;
        
    % Display
        ii = 19 ;
        clf ;
            im = imagesc(IMG(ii)) ;
            axis ij
            axis tight
            axis equal
            scat = scatter(PPP(cMASK(:),1,ii),PPP(cMASK(:),2,ii),10,sqrt(sum(UUU(cMASK(:),:,ii).^2,2))) ;
            scat.Marker = '.' ;
            
            
%% SUBSTRACTION OF THE MACRO FIELD

    % Reshape interpolated micro position to get it by quad
        PPPq = reshape(PPP,nQuads,[],2,hd.nFrames) ; % [nQuads x nPts x 2 x nFrames]
    
    % Substract the macro field
        PPPq = PPPq-PPPi ;
        
    % Restore the approximate shape
        PPPq = PPPq + cat(3,ex*Lq(1),ey*Lq(2)) ;
        
    % Compute the displacement residue (micro-macro)
        UUUq = PPPq - PPPq(:,:,:,refFrame) ;
        
    % Display
        ii = 28 ;
        qq = round(nQuads/2) ;
        ax = gobjects(0) ;
        scat = gobjects(0) ;
        clf ;
            ax(1) = mysubplot(1,2,1) ;
                scat(1) = scatter(PPPq(qq,:,1,ii),PPPq(qq,:,2,ii),400,UUUq(qq,:,1,ii)) ;
                colorbar
            ax(2) = mysubplot(1,2,2) ;
                scat(2) = scatter(PPPq(qq,:,1,ii),PPPq(qq,:,2,ii),400,UUUq(qq,:,2,ii)) ;
                colorbar
        axis(ax,'ij')
        axis(ax,'tight')
        axis(ax,'equal')
        set(scat,'Marker','.') ;
        %set(ax,'xtick',[],'ytick',[],'xcolor','none','ycolor','none')
        
%% Display the full residue map

        ii = 18 ;
        xdata = PPP(:,1,ii) ;
        ydata = PPP(:,2,ii) ;
        cdata = UUUq(:,:,1,ii) ;
        clf ;
            im = imagesc(IMG(ii)) ;
            scat = scatter(xdata(:),ydata(:),100,cdata(:)) ;
        ax = gca ;
        colorbar(ax)
        axis(ax,'ij')
        axis(ax,'tight')
        axis(ax,'equal')
        set(scat,'Marker','.') ;
        %set(ax,'xtick',[],'ytick',[],'xcolor','none','ycolor','none')
            
        
%% COMPUTE THE MACRO STRAINS

    % Derivatives Di
        % PPPi = Pq(:,1,:,:).*(1-ex).*(1-ey) ...
        %         + Pq(:,2,:,:).*(ex).*(1-ey) ...
        %         + Pq(:,3,:,:).*(ex).*(ey) ...
        %         + Pq(:,4,:,:).*(1-ex).*(ey) ;
        % Derivarives w.r.t local coordinates [ex ey]
            dx_dex = Pq(:,1,:,:).*-(1-ey) ...
                    + Pq(:,2,:,:).*(1-ey) ...
                    + Pq(:,3,:,:).*(ey) ...
                    + Pq(:,4,:,:).*-(ey) ;
            dx_dey = Pq(:,1,:,:).*-(1-ex) ...
                    + Pq(:,2,:,:).*-(ex) ...
                    + Pq(:,3,:,:).*(ex) ...
                    + Pq(:,4,:,:).*(1-ex) ;
        % Determinant of the Jacobian
            detJ = dx_dex(:,:,1,refFrame).*dx_dey(:,:,2,refFrame) - dx_dex(:,:,2,refFrame).*dx_dey(:,:,1,refFrame) ;
        % Derivatives w.r.t glocal coordinates [X1 X2]
            dx_dX1 = (1./detJ).*(dx_dey(:,:,2,refFrame).*dx_dex - dx_dex(:,:,2,refFrame).*dx_dey) ;
            dx_dX2 = (1./detJ).*(-dx_dey(:,:,1,refFrame).*dx_dex + dx_dex(:,:,1,refFrame).*dx_dey) ;
            
    % Gradient
        F11 = dx_dX1(:,:,1,:) ;
        F12 = dx_dX2(:,:,1,:) ;
        F21 = dx_dX1(:,:,2,:) ;
        F22 = dx_dX2(:,:,2,:) ;
            
    % Green-Lagrange
        L11 = 0.5*(F11.*F11 + F21.*F21 - 1) ;
        L22 = 0.5*(F12.*F12 + F22.*F22 - 1) ;
        L12 = 0.5*(F11.*F12 + F21.*F22) ;
        
    %% Display
        clf ;
            ii = 30 ;
            im = imagesc(IMG(ii)) ;
                axis ij
                axis tight
                axis equal
                xdata = PPP(:,1,ii) ;
                ydata = PPP(:,2,ii) ;
                cdata = L12(:,:,:,ii) ;
                scat = scatter(xdata(:),ydata(:),400,cdata(:)) ;
                scat.Marker = '.' ;
                colorbar
        
        
%% IDENTIFY STRAIN MODES u_ij
% u(q,x) = L_11(q)*u_11(x) + L_22(q)*u_22(x) + L_12(q)*u_12(x)
% u_ij : three modes for each frame (non linear)
    
    % Cull points not in the common mask
        uuu = UUUq(:,tMASK,:,:) ;
        l11 = L11(:,tMASK,:,:) ;
        l22 = L22(:,tMASK,:,:) ;
        l12 = L12(:,tMASK,:,:) ;
        
    % Reshape the data [nQuads nFrames nPts nDims]
        uu = permute(uuu,[1 4 2 3]) ;
        lij = permute(cat(3,l11,l22,l12),[1 4 2 3]) ;
            
    % A.u_ij = b
        fr = 30 ;
        u_ij = NaN(nCommonPts,2,3) ;
        for pp = 1:nCommonPts
            A = reshape(lij(:,fr,pp,:),[],3) ;
            b = reshape(uu(:,fr,pp,:),[],2) ;
            X = A\b ;
            u_ij(pp,:,:) = permute(X,[3 2 1]) ;
        end
        
    % Display
        clf ;
        ax = gobjects(0) ;
        scat = gobjects(0) ;
        amp = 1/20 ;
        clf ;
        for ec = 1:3
            for uc = 1:2
                ax(end+1) = mysubplot(2,3,(uc-1)*3 + ec) ;
                scat(end+1) = scatter(ax(end),...
                                ex(tMASK)'*Lq(1)+amp*u_ij(:,1,ec),...
                                ey(tMASK)'*Lq(2)+amp*u_ij(:,2,ec),...
                                400,...
                                u_ij(:,uc,ec)) ;
                colorbar(ax(end))
            end
        end
        axis(ax,'ij')
        axis(ax,'tight')
        axis(ax,'equal')
        set(scat,'Marker','.') ;
        set(ax,'xtick',[],'ytick',[],'xcolor','none','ycolor','none')
        
%% PLOT RESIDUE MAP
    ii = 30 ;
    UUUij = NaN(1,nPx*nPy,2,3) ;
    UUUij(:,tMASK,:,:) = u_ij ;
    UUUr = UUUij(:,:,:,1).*L11(:,:,:,ii) + UUUij(:,:,:,2).*L22(:,:,:,ii) + UUUij(:,:,:,3).*L12(:,:,:,ii) ;
    RRR = UUUr-UUUq(:,:,:,ii) ;
    clf ;
        im = imagesc(IMG(ii)) ;
            axis ij
            axis tight
            axis equal
            xdata = PPP(:,1,ii) ;
            ydata = PPP(:,2,ii) ;
            cdata = sum(UUUq(:,:,:,ii).^2,3) ; min(1,sum(RRR.^2,3)./sum(UUUq(:,:,:,ii).^2,3)) ;
            scat = scatter(xdata(:),ydata(:),400,cdata(:)) ;
            scat.Marker = '.' ;
            colorbar
    
        
        
