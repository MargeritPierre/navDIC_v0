%% COMPUTE CONSTANT OBJECTS

% Derivation Kernels
    N = 3 ;
    xx = (-N:N)' ;
    % Various kernel choices
        kern = [0 1 0]' ; dkern = [1 0 -1]'/2 ;
        %kern = cos(pi/2*xx/N).^2 ; dkern = -pi/N*sin(pi/2*xx/N).*cos(pi/2*xx/N) ;
        %sig = (N+1)/log(5*N) ; kern = exp(-xx.^2/sig^2) ; dkern =  -2*xx/sig^2.*exp(-xx.^2/sig^2) ;
    NORM = sum(sum(kern*kern')) ;
    Func = @(img)conv2(double(img),kern*kern','same')/NORM/imgClassRange ;
    dFunc_dx = @(img)conv2(double(img),kern*dkern','same')/NORM/imgClassRange ;
    dFunc_dy = @(img)conv2(double(img),dkern*kern','same')/NORM/imgClassRange ;

% Reference Image Processing
    % Convolution
        F = Func(img0) ;
    % Gradient of Reference image
        dF_dx = dFunc_dx(img0) ;
        dF_dy = dFunc_dy(img0) ;
        
% Projection of the Gradient in the nodal basis
    dF_da = [dF_dx(:).*MAPPING , dF_dy(:).*MAPPING] ;
    
% Hessian
    Hess = dF_da'*dF_da ;