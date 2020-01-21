%% ESTIMATE STRESS USING UNIAXIAL POWER HARDENING ELASTOPLASTIC CONSTITUTIVE LAW
% MODEL: S-Sy = K.Ep^n, SIG = E*Ee, EPS = Ee + Ep
% Plane stresses, incompressible material
% varargin : E,nu,Sy,K,n
function [S11,S22,S12] = Isotropic_Hardening(D11,D22,D12,varargin)

    % Default Parameters
        if nargin<4 ; E = 15*10e9 ; else ; E = varargin{1} ; end
        if nargin<5 ; nu = 10.4 ; else ; nu = varargin{2} ; end
        if nargin<6 ; Sy = 320e6 ; else ; Sy = varargin{3} ; end
        if nargin<7 ; K = 880e6 ; else ; K = varargin{4} ; end
        if nargin<8 ; n = 0.67 ; else ; n = varargin{5} ; end
        
    % Hard parameters
        sz = size((D11+D22+D12)) ;
        dimTime = numel(sz) ;
        nT = sz(dimTime) ;
        
    % Reshaping
        indDims = [dimTime setdiff(1:numel(sz),dimTime)] ;
        D11 = reshape(permute(D11,indDims),nT,[]) ;
        D22 = reshape(permute(D22,indDims),nT,[]) ;
        D12 = reshape(permute(D12,indDims),nT,[]) ;
        nP = size(D11,2) ;
        
    % Elastic compliance matrix SSe (in Plane Stresses)
        G = E/2/(1+nu) ;
        Es = E./(1-nu.^2) ;
        CCe = Es*[1 nu 0 ; nu 1 0 ; 0 0 0] + G*[0 0 0 ; 0 0 0 ; 0 0 1] ;
        SSe = inv(CCe) ;
        
    % Initialization
        S11 = NaN(nT,nP) ; S22 = S11 ; S12 = S11 ;
        Ep = zeros(1,nP) ; % cummulated plastic strain 
        R = Sy ; % Variable Elastic domain "radius"
        H = n*K ; % Variable Isotropic hardening modulus
            
    % Computation
        wtbr = waitbar(0,'Isotropic Hardening Stresses Computation...') ;
        for tt = 1:nT
            % Strain rate
                D = [D11(tt,:) ; D22(tt,:) ; 2*D12(tt,:)] ;
            % Equivalent strain rate (Dmean=0, incompressible)
                Deq = sqrt(2/3*(D(1,:).^2 + D(2,:).^2 + 1/2*D(3,:).^2)) ;
            % First Elastic step
                if tt==1
                    S11(tt,:) = CCe(1,:)*D ;
                    S22(tt,:) = CCe(2,:)*D ;
                    S12(tt,:) = CCe(3,:)*D ;
                    continue
                end
            % STRESS INCREMENT
                % Mean stress
                    Smean = 1/3.*(S11(tt-1,:)+S22(tt-1,:)) ;
                % Deviatoric stresses
                    s11 = S11(tt-1,:)-Smean ;
                    s22 = S22(tt-1,:)-Smean ;
                    s33 = -Smean ;
                    s12 = S12(tt-1,:) ;
                % Equivalent VM stress
                    Seq = sqrt(3/2*(s11.^2 + s22.^2 + s33.^2 + 2*s12.^2)) ;
                % Yield surface Normal
                    N11 = 3/2*s11./Seq ;
                    N22 = 3/2*s22./Seq ;
                    N12 = 3/2*s12./Seq ;
                % Plastic tangent compliance (Sijkl=NijNkl)
                    SSp1111 = N11.^2./H ;
                    SSp1112 = N11.*N12./H ;
                    SSp1122 = N11.*N22./H ;
                    SSp1212 = N12.^2./H ;
                    SSp1222 = N12.*N22./H ;
                    SSp2222 = N22.*N22./H ;
                % In voigt Notation...
                    SSp11 = SSp1111 ;
                    SSp12 = SSp1122 ;
                    SSp16 = 2*SSp1112 ;
                    SSp22 = SSp2222 ;
                    SSp26 = 2*SSp1222 ;
                    SSp66 = 4*SSp1212 ;
                % As a matrix
                    SSp = reshape([SSp11 ; SSp12 ; SSp16 ; SSp12 ; SSp22 ; SSp16 ; SSp16 ; SSp26 ; SSp66],3,3,nP) ;
                % Yield criterion
                    isYielding = Seq>R ;
                    SSp(:,:,~isYielding) = 0 ;
                % Tangent Compliance
                    SSt = SSe + SSp ;
                % Tangent Stiffness
                    CCt = zeros(size(SSt)) ;
                    for pp = 1:nP
                        CCt(:,:,pp) = inv(SSt(:,:,pp)) ;
                    end
                % Stress increment
                    dS = reshape(sum(CCt.*reshape(D,1,3,nP),2),3,nP) ;
                    S11(tt,:) = S11(tt-1,:) + dS(1,:) ;
                    S22(tt,:) = S22(tt-1,:) + dS(2,:) ;
                    S12(tt,:) = S12(tt-1,:) + dS(3,:) ;
            % YIELD SURFACE UPDATE
                if all(~isYielding) ; continue ; end
                % Flow rule
                    hasYielded = N11.*dS(1,:) + N22.*dS(2,:) + 2*N12.*dS(3,:) > 0 ;
                    hasYielded(isYielding) = true ;
                % Cummulated plastic strain
                    Ep(hasYielded) = Ep(hasYielded) + abs(Deq(hasYielded)) ;
                % Elastic domain size
                    R = K.*Ep.^n + Sy ;
                % Isotropic Hardening modulus dR/dp
                    H = n*K.*Ep.^(n-1) ;
            % Waitbar
                wtbr = waitbar(tt/nT,wtbr) ;
        end
        delete(wtbr) ;
        
    % Re-Reshape
        [~,indRev] = sort(indDims) ;
        S11 = permute(reshape(S11,sz(indDims)),indRev) ;
        S22 = permute(reshape(S22,sz(indDims)),indRev) ;
        S12 = permute(reshape(S12,sz(indDims)),indRev) ;

return


%% TEST THE MODEL

    E = 150e9 ;
    nu = 0.5 ;
    Sy = 320e6 ;
    K = 880e6 ;
    n = .67 ; %.67 ;
    tol = abs(1000*eps*Sy) ;
    
    ax = gca ;
    EPS = linspace(0,ax.XLim(2),10000) ;
    de = 0.00005 ; EPS = [0:de:0.2 0.2:-de:0.198 0.198:de:0.3 0.3:-de:0.298 0.298:de:0.35] ;
    D = [0 diff(EPS)] ;
    D11 = D ; D22 = -0.5*D ; D12 = D.*0 ;
    
    [S11,S22,S12] = Isotropic_Hardening(D11,D22,D12,E,nu,Sy,K,n) ;
    
    Smean = 1/3*(S11+S22) ;
    SIG = sqrt(3/2*((S11-Smean).^2 + (S22-Smean).^2 + Smean.^2 + 2*S12.^2)) ;

    tag = 'Power Hardening' ;
    delete(findobj(gca,'tag',tag)) ;
    plot(EPS,SIG*1e-6,'-.k','linewidth',1,'tag',tag,'DisplayName','PH') ;
    
    
