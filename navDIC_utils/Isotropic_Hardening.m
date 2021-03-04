%% ESTIMATE STRESS USING UNIAXIAL POWER HARDENING ELASTOPLASTIC CONSTITUTIVE LAW
% MODEL: S-Sy = K.Ep^n, SIG = E*Ee, EPS = Ee + Ep
% Plane stresses, incompressible material
% varargin : E,nu,Sy,K,n,tol,maxIt
function [S11,S22,S33,S12,Dp11,Dp22,Dp33,Dp12,R] = Isotropic_Hardening(D11,D22,D12,varargin)

    % Default Parameters
        if nargin<4 ; E = 15*10e9 ; else ; E = varargin{1} ; end
        if nargin<5 ; nu = 0.5 ; else ; nu = varargin{2} ; end
        if nargin<6 ; Sy = 320e6 ; else ; Sy = varargin{3} ; end
        if nargin<7 ; K = 880e6 ; else ; K = varargin{4} ; end
        if nargin<8 ; n = 0.67 ; else ; n = varargin{5} ; end
        if nargin<9 ; tol = Sy*1e-6 ; else ; tol = varargin{6} ; end
        if nargin<10 ; maxIt = 100 ; else ; maxIt = varargin{7} ; end
        incompressible = true ;
        planeStresses = true ;
        
    % Hard parameters
        sz = size((D11+D22+D12)) ;
        dimTime = find(sz~=1,1,'last') ;
        nT = sz(dimTime) ;
        
    % Reshaping
        indDims = [setdiff(1:numel(sz),dimTime) dimTime] ;
        D11 = reshape(permute(D11,[numel(sz)+1 indDims]),1,[],nT) ;
        D22 = reshape(permute(D22,[numel(sz)+1 indDims]),1,[],nT) ;
        D12 = reshape(permute(D12,[numel(sz)+1 indDims]),1,[],nT) ;
        if incompressible % INCOMPRESSIBLE MATERIAL !
            D = [D11 ; D22 ; -D11-D22 ; D12] ;
        else
            D = [D11 ; D22 ; 2*D12] ;
        end
        nP = size(D,2) ;
        
    % Lame parameters
        lambda = E*nu/(1+nu)/(1-2*nu) ;
        mu = E/2/(1+nu) ;
        
    % Initialization
        S = D*NaN ;
        Dp = D*NaN ; % Plastic Strain Rate
        R = D(1,:,:)*NaN ; % Variable Elastic domain "radius"
            
    % Computation
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
                Str = S(:,:,tt) + 2*mu*D(:,:,tt) ;
            % STRESS DEVIATOR
                Smean = 1/3*sum(Str(1:3,:),1) ;
                Sdev = Str - [1;1;1;0].*Smean ;
            % EQUIVALENT STRESS
                Seq = sqrt(3/2*sum([1;1;1;2].*(Sdev.^2),1)) ;
            % ELASTICITY DOMAIN NORMAL
                N = 3/2*Sdev./(Seq+eps) ; % 'eps' is here to avoid NaN when Seq vanishes
            % PLASTIC CORRECTION
                P = sum(sqrt(2/3*sum([1;1;1;2].*(Dp(:,:,1:tt).^2),1)),3) ;
                p = P ;
                for it = 1:maxIt
                    % Yield Criterion
                        f = Seq-3*mu*(p-P)-R(:,:,tt) ;
                        yield = f>tol ;
                        if all(~yield) ; break ; end
                    % Update the cummulated plastic strain increment
                        dR_dp = n*K*(p(yield)+eps).^(n-1) ; % eps is here to prevent division by zero
                        df_dp = - (3*mu + dR_dp) ;
                        p(yield) = p(yield) - f(yield)./df_dp ;
                    % Update the elastic domain
                        R(:,:,tt) = Sy + K*p.^n ;
                end
            % PLASTIC FLOW
                dp = p-P ;
                Dp(:,:,tt) = dp.*N ;
            % STRESS COMPUTATION
                Ee = sum(D(:,:,1:tt)-Dp(:,:,1:tt),3) ;
                S(:,:,tt) = 2*mu*(Ee - [1;1;1;0].*Ee(3,:)) ;
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
        S33 = reshp(S(3,:,:)) ;
        S12 = reshp(S(4,:,:)) ;
        R = reshp(R) ;
        Dp11 = reshp(Dp(1,:,:)) ;
        Dp22 = reshp(Dp(2,:,:)) ;
        Dp33 = reshp(Dp(3,:,:)) ;
        Dp12 = reshp(Dp(4,:,:)) ;

return


%% TEST THE MODEL

    E = 200e9 ;
    nu = .5 ;
    Sy = 50*10*1e6 ;
    K = 167*10*1e6 ;
    n = 1 ;
    tol = abs(1000*eps*Sy) ;
    
    ax = gca ;
    EPS = linspace(0,ax.XLim(2),1000) ;
    %de = 0.0005 ; EPS = [0:de:0.2 0.2:-de:0.1 0.1:de:0.3 0.3:-de:0.298 0.298:de:0.35] ;
    de = 0.0001 ; EPS = [0:de:0.335 0.335:-de:0.333] ;
    D = [0 diff(EPS)]' ;
    D11 = D ; D22 = -nu*D ; D12 = D.*0 ;
    
    [S11,S22,S33,S12,Dp11,Dp22,Dp33,Dp12,R] = Isotropic_Hardening(D11,D22,D12,E,nu,Sy,K,n) ;
    
    Smean = 1/3*(S11+S22) ;
    SIG = sqrt(3/2*((S11-Smean).^2 + (S22-Smean).^2 + Smean.^2 + 2*S12.^2)) ;

    tag = 'Power Hardening' ;
    delete(findobj(gca,'tag',tag)) ;
    plot(EPS,S11*1e-6,'-k','linewidth',2,'tag',tag,'DisplayName','PH') ;
    %plot(EPS,R*1e-6,'--r','linewidth',1,'tag',tag,'DisplayName','R') ;
    %plot(EPS,(Sy+K.*(EPS(:)-SIG(:)/E).^n)*1e-6,':b','linewidth',1,'tag',tag,'DisplayName','Rth') ;
    ax.XLim(1) = 0 ; ax.YLim(1) = 0 ; 
    
    
    
    
    
    
    
    
    
return
    
