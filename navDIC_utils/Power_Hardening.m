%% ESTIMATE STRESS USING UNIAXIAL POWER HARDENING ELASTOPLASTIC CONSTITUTIVE LAW
% MODEL: S-Sy = K.Ep^n, SIG = E*Ee, EPS = Ee + Ep
function [SIG,Es,Nus] = Power_Hardening(EPS,E,Sy,K,n,Nu)

    % Default Parameters
        if nargin==1
            E = 15*10e9 ;
            Sy = 320e6 ;
            K = 880e6 ;
            n = 0.67 ;
            Nu = 0.4;
        end
        

    % Initialize Stresses
        SIG = Sy*ones(size(EPS)) ; 
    % Elastic regime
        Ey = Sy/E ;
        isElastic = abs(EPS)<Ey ;
        SIG(isElastic) = E*EPS(isElastic) ;
        SIG(~isElastic) = (Sy+K.*(abs(EPS(~isElastic))-Ey).^n).*sign(EPS(~isElastic)) ;
        
    % Modulus and Poisson Ratio
        if nargout>1
            Es = SIG./EPS ;
            Es(SIG==0 & EPS==0) = E ;
            Nus = 0.5 + Es/E*(Nu-0.5) ; 
        end

return


%% POWER HARDING MODEL

    E = 150e9 ;
    Sy = 320e6 ;
    K = 880e6 ;
    n = .67 ;
    Nu = 0.33;
    tol = abs(1000*eps*Sy) ;
    
    ax = gca ;
    EPS = linspace(0,ax.XLim(2),10000) ;
    
    [SIG,Es,Nus] = Power_Hardening(EPS,tol,E,Sy,K,n,Nu) ;

    tag = 'Power Hardening' ;
    delete(findobj(gca,'tag',tag)) ;
    plot(EPS,SIG*1e-6,'-.k','linewidth',1,'tag',tag,'DisplayName','PH') ;
    
    
