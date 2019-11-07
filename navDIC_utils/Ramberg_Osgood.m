%% ESTIMATE STRESS USING UNIAXIAL RAMBERG-OSGOOD CONSTITUTIVE LAW
% MODEL: EPS = SIG/E.*(1+alpha*(SIG/Sy).^(n-1))
function [SIG,Es,Nus] = Ramberg_Osgood(EPS,tol,E,Sy,n,alpha,Nu)

    % Default Parameters
        if nargin==1
            E = 18*10e9 ;
            Sy = 300e6 ;
            n = 91/10 ;
            alpha = E/Sy*0.002 ; 
            Nu = 0.33;
            tol = abs(1000*eps*Sy) ;
        end

    % Initialize Stresses
        SIG = Sy*ones(size(EPS)) ; 
        dSIG = Inf ;
    % Newton-Raphson procedure
        it = 0 ;
        while max(abs(dSIG(:)))>tol
            % Newton-Raphson increment
                dEPS = SIG/E.*(1+alpha*(SIG/Sy).^(n-1)) - EPS ;
                dEPS_dSIG = 1/E.*(1+alpha*(SIG/Sy).^(n-1)) + SIG/E.*(alpha/Sy*(n-1).*(SIG/Sy).^(n-2)) ;
                dSIG = -dEPS./dEPS_dSIG ;
            % Before yield
                SIG = SIG + dSIG ;
            % Convergence criterion
                it = it + 1 ;
        end
        
    % Modulus and Poisson Ratio
        if nargout>1
            Es = E./(1+alpha*(SIG/Sy).^(n-1));
            Nus = 0.5+Es/E*(Nu-0.5);
        end

end


%% RAMBERG-OSGOOD MODEL

%     E = 18*10e9 ;
%     Sy = 300e6 ;
%     n = 91/10 ;
%     alpha = E/Sy*0.002 ;
% 
%     SIG = linspace(0,600e6,1000) ;
%     EPS = SIG/E.*(1+alpha*(SIG/Sy).^(n-1)) ;
%     Ep = SIG/E.*(1+alpha*(SIG/Sy).^(n-1))
% 
%     tag = 'jnfeubue' ;
%     delete(findobj(gca,'tag',tag)) ;
%     plot(EPS,SIG*1e-6,'k','tag',tag) ;
    
    
