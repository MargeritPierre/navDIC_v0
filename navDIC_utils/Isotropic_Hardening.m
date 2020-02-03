%% ESTIMATE STRESS USING UNIAXIAL POWER HARDENING ELASTOPLASTIC CONSTITUTIVE LAW
% MODEL: S-Sy = K.Ep^n, SIG = E*Ee, EPS = Ee + Ep
% Plane stresses, incompressible material
% varargin : E,nu,Sy,K,n
function [S11,S22,S12,Dp11,Dp22,Dp12,R] = Isotropic_Hardening(D11,D22,D12,varargin)

    % Default Parameters
        if nargin<4 ; E = 15*10e9 ; else ; E = varargin{1} ; end
        if nargin<5 ; nu = 0.4 ; else ; nu = varargin{2} ; end
        if nargin<6 ; Sy = 320e6 ; else ; Sy = varargin{3} ; end
        if nargin<7 ; K = 880e6 ; else ; K = varargin{4} ; end
        if nargin<8 ; n = 0.67 ; else ; n = varargin{5} ; end
        if nargin<9 ; tol = Sy*1e-6 ; else ; tol = varargin{6} ; end
        if nargin<10 ; maxIt = 100 ; else ; maxIt = varargin{7} ; end
        
    % Hard parameters
        sz = size((D11+D22+D12)) ;
        dimTime = find(sz~=1,1,'last') ;
        nT = sz(dimTime) ;
        
    % Reshaping
        indDims = [setdiff(1:numel(sz),dimTime) dimTime] ;
        D11 = reshape(permute(D11,[numel(sz)+1 indDims]),1,[],nT) ;
        D22 = reshape(permute(D22,[numel(sz)+1 indDims]),1,[],nT) ;
        D12 = reshape(permute(D12,[numel(sz)+1 indDims]),1,[],nT) ;
        D = [D11 ; D22 ; 2*D12] ;
        nP = size(D,2) ;
        
    % Deviatoric projection matrix
        Pd = 1/3*[2 -1 0 ; -1 2 0 ; 0 0 6] ;
        
    % Elastic compliance matrix SSe (in Plane Stresses)
        G = E/2/(1+nu)./(1-nu.^2) ;
        Es = E./(1-nu.^2) ;
        C = Es*[1 nu 0 ; nu 1 0 ; 0 0 0] + G*[0 0 0 ; 0 0 0 ; 0 0 1] ;
        
    % Initialization
        S = D*NaN ;
        Dp = D*NaN ; % Plastic Strain Rate
        R = D(1,:,:)*NaN ; % Variable Elastic domain "radius"
            
    % Computation
    % see 'Return Mapping for Plane Stresses' by Simo&Hugues's book, p.130
        tik = tic ;
        wtbr = waitbar(0,'Isotropic Hardening Stresses Computation...') ;
        for tt = 1:nT
            % INITIALIZATION
                if tt==1 % First step
                    S(:,:,tt) = 0 ;
                    Dp(:,:,tt) = 0 ;
                    R(:,:,tt) = Sy ;
                else
                    S(:,:,tt) = S(:,:,tt-1) ;
                    Dp(:,:,tt) = 0 ;
                    R(:,:,tt) = R(:,:,tt-1) ;
                end
            % ELASTIC PREDICTOR
                S(:,:,tt) = S(:,:,tt) + C*D(:,:,tt) ;
            for it = 1:maxIt
            % PLASTIC CORRECTION
                % Equivalent Stress
                    Seq = sqrt(3/2*sum(S(:,:,tt).*(Pd*S(:,:,tt)),1)) ;
                % Yield Criterion
                    f = Seq-R(:,:,tt) ;
                    yield = f>tol ;
                    if all(~yield) ; break ; end
                % Yield Surface Normal
                    N = 3/2*(Pd*S(:,yield,tt))./Seq(yield) ;
                % Flow rule
                    lmbda = f(yield)./sum(N.*(C*N),1) ;
                % Plastic Flow
                    dEp = lmbda.*(N) ;
                    Dp(:,yield,tt) = Dp(:,yield,tt) + dEp ;
                % Plastic Correction
                    S(:,yield,tt) = S(:,yield,tt) - C*dEp ;
            % YIELD SURFACE UPDATE
                % Cummulated eq. plastic strains
                    Ep = sum(sqrt(2/3*(Dp(1,yield,1:tt).^2 + Dp(2,yield,1:tt).^2 + 1/2*Dp(3,yield,1:tt).^2)),3) ;
                % Elastic domain radius
                    R(:,yield,tt) = Sy + K.*Ep.^n ;
            end
            % Waitbar
                if toc(tik)>0.1
                    wtbr = waitbar(tt/nT,wtbr) ;
                    tik = tic ;
                end
        end
        delete(wtbr) ;

        
    % Re-Reshape
        [~,indRev] = sort(indDims) ;
        reshp = @(data)permute(reshape(data,sz(indDims)),indRev) ;
        S11 = reshp(S(1,:,:)) ;
        S22 = reshp(S(2,:,:)) ;
        S12 = reshp(S(3,:,:)) ;
        R = reshp(R) ;
        Dp11 = reshp(Dp(1,:,:)) ;
        Dp22 = reshp(Dp(2,:,:)) ;
        Dp12 = 0.5*reshp(Dp(3,:,:)) ;

return


%% TEST THE MODEL

    E = 150e9 ;
    nu = .5 ;
    Sy = 320e6 ;
    K = 880e6 ;
    n = 1 ; %.67 ;
    tol = abs(1000*eps*Sy) ;
    
    ax = gca ;
    EPS = linspace(0,ax.XLim(2),10000) ;
    de = 0.0001 ; EPS = [0:de:0.2 0.2:-de:0.1 0.1:de:0.3 0.3:-de:0.298 0.298:de:0.35] ;
    D = [0 diff(EPS)]' ;
    D11 = D ; D22 = -nu*D ; D12 = D.*0 ;
    
    [S11,S22,S12,Dp11,Dp22,Dp12,R] = Isotropic_Hardening(D11,D22,D12,E,nu,Sy,K,n) ;
    
    Smean = 1/3*(S11+S22) ;
    SIG = sqrt(3/2*((S11-Smean).^2 + (S22-Smean).^2 + Smean.^2 + 2*S12.^2)) ;

    tag = 'Power Hardening' ;
    delete(findobj(gca,'tag',tag)) ;
    plot(EPS,S11*1e-6,'.-.k','linewidth',1,'tag',tag,'DisplayName','PH') ;
    plot(EPS,R*1e-6,'--r','linewidth',1,'tag',tag,'DisplayName','R') ;
    plot(EPS,(Sy+K.*(EPS-SIG/E).^n)*1e-6,':b','linewidth',1,'tag',tag,'DisplayName','Rth') ;
    
    
    
    
    
    
    
    
    
return







%% THIRD ATTEMPT
%     tik = tic ;
%     wtbr = waitbar(0,'Isotropic Hardening Stresses Computation...') ;
%     for tt = 1:nT
%         % INITIALIZATION
%             if tt==1 % First step
%                 S(:,:,tt) = 0 ;
%                 Dp(:,:,tt) = 0 ;
%                 R(:,:,tt) = Sy ;
%             else
%                 S(:,:,tt) = S(:,:,tt-1) ;
%                 Dp(:,:,tt) = 0 ;
%                 R(:,:,tt) = R(:,:,tt-1) ;
%             end
%         % ELASTIC PREDICTOR (Trial Stresses)
%             St = C*(sum(D(:,:,1:tt)-Dp(:,:,1:tt),3)) ;
%         % YIELD CRITERION CHECK
%             % Equivalent Stress
%                 Seq = sqrt(3/2*sum(St.*(Pd*St),1)) ;
%             % Yield Criterion
%                 f = Seq-R(:,:,tt) ;
%                 yield = f>tol ;
%                 if all(~yield) ; break ; end
%         % PLATICITY MULTIPLIER COMPUTATION
%             if any(yield)
%                 % Current cummulative plastic strain
%                     Ep = sum(sqrt(2/3*Dp(:,yield,1:tt)).^2,3) ;
%                 % Functions
%                     fbar2 = @(gamma,yield) 1/6*(St(1,yield)+St(2,yield)).^2./(1+E/3/(1-nu).*gamma).^2 + 1/2*((St(1,yield)-St(2,yield)).^2 + 2*St(3,yield).^2)./(1+2*G.*gamma).^2 ;
%                     Rbar2 = @(gamma,yield) (Sy+K.*(Ep(yield)+gamma.*sqrt(fbar2(gamma,yield))).^n).^2 ;
%                     f2 = @(gamma,yield) fbar2(gamma,yield)- Rbar2(gamma,yield) ;
%                 % Minimize (constitency)
%                     for it = 1:maxIt
%                         % Yield Surface Normal
%                             N = 3/2*(Pd*S(:,yield,tt))./Seq(yield) ;
%                         % Flow rule
%                             lmbda = f./sum(N.*(C*N),1) ;
%                         % Plastic Flow
%                             dEp = lmbda.*N ;
%                             Dp(:,yield,tt) = Dp(:,yield,tt) + dEp ;
%                         % Plastic Correction
%                             S(:,yield,tt) = S(:,yield,tt) - C*dEp ;
%                     % YIELD SURFACE UPDATE
%                         % Cummulated eq. plastic strains
%                             Ep = sum(sqrt(2/3*(Dp(1,yield,1:tt).^2 + Dp(2,yield,1:tt).^2 + 1/2*Dp(3,yield,1:tt).^2)),3) ;
%                         % Elastic domain radius
%                             R(:,yield,tt) = Sy + K.*Ep.^n ;
%                     end
%             end
%         % Waitbar
%             if toc(tik)>0.1
%                 wtbr = waitbar(tt/nT,wtbr) ;
%                 tik = tic ;
%             end
%     end
%     delete(wtbr) ;



%% SECOND ATTEMPT
%         tik = tic ;
%         wtbr = waitbar(0,'Isotropic Hardening Stresses Computation...') ;
%         for tt = 1:nT
%             % INITIALIZATION
%                 if tt==1 % First step
%                     S(:,:,tt) = 0 ;
%                     Dp(:,:,tt) = 0 ;
%                     R(:,:,tt) = Sy ;
%                 else
%                     S(:,:,tt) = S(:,:,tt-1) ;
%                     Dp(:,:,tt) = 0 ;
%                     R(:,:,tt) = R(:,:,tt-1) ;
%                 end
%             % ELASTIC PREDICTOR
%                 S(:,:,tt) = S(:,:,tt) + C*D(:,:,tt) ;
%             for it = 1:maxIt
%             % PLASTIC CORRECTION
%                 % Equivalent Stress
%                     Seq = sqrt(3/2*sum(S(:,:,tt).*(Pd*S(:,:,tt)),1)) ;
%                 % Yield Criterion
%                     f = Seq-R(:,:,tt) ;
%                     yield = f>tol ;
%                     if all(~yield) ; break ; end
%                 % Yield Surface Normal
%                     N = 3/2*(Pd*S(:,yield,tt))./Seq(yield) ;
%                 % Flow rule
%                     lmbda = f./sum(N.*(C*N),1) ;
%                 % Plastic Flow
%                     dEp = lmbda.*N ;
%                     Dp(:,yield,tt) = Dp(:,yield,tt) + dEp ;
%                 % Plastic Correction
%                     S(:,yield,tt) = S(:,yield,tt) - C*dEp ;
%             % YIELD SURFACE UPDATE
%                 % Cummulated eq. plastic strains
%                     Ep = sum(sqrt(2/3*(Dp(1,yield,1:tt).^2 + Dp(2,yield,1:tt).^2 + 1/2*Dp(3,yield,1:tt).^2)),3) ;
%                 % Elastic domain radius
%                     R(:,yield,tt) = Sy + K.*Ep.^n ;
%             end
%             % Waitbar
%                 if toc(tik)>0.1
%                     wtbr = waitbar(tt/nT,wtbr) ;
%                     tik = tic ;
%                 end
%         end
%         delete(wtbr) ;




%% OLD IMPLEMENTATION (TANGENT STIFFNESS)
            
    % Computation
%         wtbr = waitbar(0,'Isotropic Hardening Stresses Computation...') ;
%         for tt = 1:nT
%             % Strain rate
%                 D = [D11(tt,:) ; D22(tt,:) ; 2*D12(tt,:)] ;
%                 if D11(tt+2,:)<0
%                     D(1)
%                 end
%             % Equivalent strain rate (Dmean=0, incompressible)
%                 Deq = sqrt(2/3*(D(1,:).^2 + D(2,:).^2 + 1/2*D(3,:).^2)) ;
%             % First Elastic step
%                 if tt==1
%                     S11(tt,:) = CCe(1,:)*D ;
%                     S22(tt,:) = CCe(2,:)*D ;
%                     S12(tt,:) = CCe(3,:)*D ;
%                     continue
%                 end
%             % STRESS INCREMENT
%                 % Mean stress
%                     Smean = 1/3.*(S11(tt-1,:)+S22(tt-1,:)) ;
%                 % Deviatoric stresses
%                     s11 = S11(tt-1,:)-Smean ;
%                     s22 = S22(tt-1,:)-Smean ;
%                     s33 = -Smean ;
%                     s12 = S12(tt-1,:) ;
%                 % Equivalent VM stress
%                     Seq = sqrt(3/2*(s11.^2 + s22.^2 + s33.^2 + 2*s12.^2)) ;
%                 % Yield surface Normal
%                     N11 = 3/2*s11./Seq ;
%                     N22 = 3/2*s22./Seq ;
%                     N12 = 3/2*s12./Seq ;
%                 % Plastic tangent compliance (Sijkl=NijNkl)
%                     SSp1111 = N11.^2./H ;
%                     SSp1112 = N11.*N12./H ;
%                     SSp1122 = N11.*N22./H ;
%                     SSp1212 = N12.^2./H ;
%                     SSp1222 = N12.*N22./H ;
%                     SSp2222 = N22.*N22./H ;
%                 % In voigt Notation...
%                     SSp11 = SSp1111 ;
%                     SSp12 = SSp1122 ;
%                     SSp16 = 2*SSp1112 ;
%                     SSp22 = SSp2222 ;
%                     SSp26 = 2*SSp1222 ;
%                     SSp66 = 4*SSp1212 ;
%                 % As a matrix
%                     SSp = reshape([SSp11 ; SSp12 ; SSp16 ; SSp12 ; SSp22 ; SSp16 ; SSp16 ; SSp26 ; SSp66],3,3,nP) ;
%                 % Yield criterion
%                     isYielding = Seq>R ;
%                     SSp(:,:,~isYielding) = 0 ;
%                 % Tangent Compliance
%                     SSt = SSe + SSp ;
%                 % Tangent Stiffness
%                     CCt = zeros(size(SSt)) ;
%                     for pp = 1:nP
%                         CCt(:,:,pp) = inv(SSt(:,:,pp)) ;
%                     end
%                 % Stress increment
%                     dS = reshape(sum(CCt.*reshape(D,1,3,nP),2),3,nP) ;
%                     S11(tt,:) = S11(tt-1,:) + dS(1,:) ;
%                     S22(tt,:) = S22(tt-1,:) + dS(2,:) ;
%                     S12(tt,:) = S12(tt-1,:) + dS(3,:) ;
%             % YIELD SURFACE UPDATE
%                 if all(~isYielding) ; continue ; end
%                 % Flow rule
%                     hasYielded = N11.*dS(1,:) + N22.*dS(2,:) + 2*N12.*dS(3,:) > 0 ;
%                     hasYielded(isYielding) = true ;
%                 % Cummulated plastic strain
%                     Ep(hasYielded) = Ep(hasYielded) + abs(Deq(hasYielded)) ;
%                 % Elastic domain size
%                     R = K.*Ep.^n + Sy ;
%                 % Isotropic Hardening modulus dR/dp
%                     H = n*K.*Ep.^(n-1) ;
%             % Waitbar
%                 wtbr = waitbar(tt/nT,wtbr) ;
%         end
    
