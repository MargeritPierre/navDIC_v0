%% IMAGE SMOOTHING AND DERIVATION KERNELS
    
    % Kenel relative position
        xx = (-sizeImageKernel:sizeImageKernel)' ;
    % Various kernel choices
        switch kernelModel 
            case 'finiteDiff' 
                kern = [0 1 0]' ; 
                dkern = [1 0 -1]'/2 ;
            case 'gaussian' 
                sig = (sizeImageKernel+1)/log(5*sizeImageKernel) ; % optimized variance
                kern = exp(-xx.^2/sig^2) ; 
                dkern =  -2*xx/sig^2.*exp(-xx.^2/sig^2) ;
            case 'cos2'
                kern = cos(pi/2*xx/sizeImageKernel).^2 ; 
                dkern = -pi/sizeImageKernel*sin(pi/2*xx/sizeImageKernel).*cos(pi/2*xx/sizeImageKernel) ;
        end
    % 2D Kernels
        KERN = kern*kern' ;
        kern = [0 1 0]' ; 
        dkern = [1 0 -1]'/2 ;
        dKERN_dx = kern*dkern' ;
        dKERN_dy = dkern*kern' ;
    % Declare smoothing and derivation functions
        NORM = sum(sum(KERN)) ;
        Smooth = @(img)conv2(double(img),kern*kern','same')/NORM/imgClassRange ;
        dI_dx = @(img)conv2(img,dKERN_dx,'same') ;
        dI_dy = @(img)conv2(img,dKERN_dy,'same') ;
    % Transfer function to convert to sparse data (x10 speed for big DIC meshes)
        spdiag = @(data) sparse(1:numel(data),1:numel(data),data(:)) ;
    % Jacobian corresponding to the degrees of freedom
        dI_da = @(img)[...
                    spdiag(dI_dx(img))*MAPPING ...
                    spdiag(dI_dy(img))*MAPPING ...
                    ] ;


