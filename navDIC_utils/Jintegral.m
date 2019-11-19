%% J-INTEGRAL ESTIMATION ON DIC DATA


%% IMPORT THE SEED AND COMPUTE STRESSES AND STRAIN ENERGY
    seedNumber = 2 ;
    
    % Import hd struct
        global hd
    
    % Re-compute the seed Data Fields on points
        hd.Seeds(seedNumber).DataOnNodes = true ;
    % Compute equivalent stresses
        EPS = hd.Seeds(seedNumber).DataFields.Leq ;
        [SIG,Es,Nus] = Ramberg_Osgood(EPS) ;
        hd.Seeds(seedNumber).DataFields.S11 = Es./(1-Nus.^2).*(hd.Seeds(seedNumber).DataFields.L11 + Nus.*hd.Seeds(seedNumber).DataFields.L22) ;
        hd.Seeds(seedNumber).DataFields.S22 = Es./(1-Nus.^2).*(hd.Seeds(seedNumber).DataFields.L22 + Nus.*hd.Seeds(seedNumber).DataFields.L11);
        hd.Seeds(seedNumber).DataFields.S12 = Es./(1-Nus.^2).*hd.Seeds(seedNumber).DataFields.L12;
        hd.Seeds(seedNumber).DataFields.Seq = SIG ;
    % Compute strain energy
        hd.Seeds(seedNumber).DataFields.W = SIG.*EPS ;

%% INITIALIZE THE FIGURE
    % Init the figure
        clf ;
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
        
%% SET THE INTEGRAL CONTOUR INTERACTIVELY AND COMPUTE THE INTEGRAL
    
    crackPos = [0.5 0.82].*[nJ nI] ; % Crack tip position
    crackVec = [0 -1] ; % Crack direction
    zoneWidth = 0.8*nJ ;
    zoneHeight = 0.68*nI ;
    crackWidth = 0.05*nJ ;
    rulerLength = 0.007 ;
    nPts = 1000 ;
    
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


% COMPUTE THE J-INTEGRAL
    
    % Units
        pixelsByMeters = sqrt(sum(diff(ruler.Position,1,1).^2))/rulerLength ;
        
    % Integration points
        dL = sqrt(sum(diff(contour.Position,1,1).^2,2)) ;
        L = cumsum([0 ; dL]) ;
        s = linspace(0,L(end),nPts) ; 
        Pts = interp1(L,contour.Position,s) ;
        PtsMid = (Pts(1:end-1,:)+Pts(2:end,:))/2 ;
        dx = diff(Pts,1,1)/pixelsByMeters ;
        dl = sqrt(sum(dx.^2,2)) ;
        
    % Outgoing Normals
        normals = dx./dl*[0 1 ; -1 0]' ;
        switched = inpolygon(PtsMid(:,1)+normals(:,1),PtsMid(:,2)+normals(:,2),contour.Position(:,1),contour.Position(:,2)) ;
        normals(switched,:) = -normals(switched,:) ;
        
    % Interpolation matrix
        refPts = hd.Seeds(seedNumber).MovingPoints(:,:,1) ;
        T = hd.Seeds(seedNumber).interpMat(PtsMid,0,refPts) ;
        
    % Interpolated Data
        W = T*squeeze(hd.Seeds(seedNumber).DataFields.W) ;
        du1_dx1 = T*squeeze(hd.Seeds(seedNumber).DataFields.F11-1) ;
        du2_dx1 = T*squeeze(hd.Seeds(seedNumber).DataFields.F21) ;
        du1_dx2 = T*squeeze(hd.Seeds(seedNumber).DataFields.F12) ;
        du2_dx2 = T*squeeze(hd.Seeds(seedNumber).DataFields.F22-1) ;
        S11 = T*squeeze(hd.Seeds(seedNumber).DataFields.S11) ;
        S22 = T*squeeze(hd.Seeds(seedNumber).DataFields.S22) ;
        S12 = T*squeeze(hd.Seeds(seedNumber).DataFields.S12) ;
        t1 = normals(:,1).*S11 + normals(:,2).*S12 ;
        t2 = normals(:,1).*S12 + normals(:,2).*S22 ;
        
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
        
%% PLOT THE RESULT

    pl = plot(J/1000) ;
    contour.Color = pl.Color ;

    xlabel 'Image'
    ylabel 'J (kJ/m$^2$)'

    box on
    grid on
    
    
%% ADD UNLOADING FIELDS

    hd.Seeds(seedNumber).DataFields.Deq = cat(3,...
                                                hd.Seeds(seedNumber).DataFields.Leq(:,:,1)*0,...
                                                diff(hd.Seeds(seedNumber).DataFields.Leq,1,3)) ;

    hd.Seeds(seedNumber).DataFields.Dmax = hd.Seeds(seedNumber).DataFields.Leq-max(hd.Seeds(seedNumber).DataFields.Leq,[],3) ;
    
    hd.Seeds(seedNumber).DataFields.logW = log10(hd.Seeds(seedNumber).DataFields.W) ;
 
    
    
    
    
    
