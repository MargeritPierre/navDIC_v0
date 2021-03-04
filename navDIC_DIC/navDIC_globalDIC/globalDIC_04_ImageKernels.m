%% IMAGE SMOOTHING AND DERIVATION KERNELS
    
    % Kenel relative position
        xx = (-sizeImageKernel:sizeImageKernel)' ;
        
    % Various kernel choices
        switch kernelModel 
            case 'finiteDiff' 
                kern = [0 1 0]' ; 
                %dkern = [1 0 -1]'/2 ;
            case 'gaussian' 
                sig = (sizeImageKernel+1)/log(5*sizeImageKernel) ; % optimized variance
                kern = exp(-xx.^2/sig^2) ; 
                %dkern =  -2*xx/sig^2.*exp(-xx.^2/sig^2) ;
            case 'cos2'
                kern = cos(pi/2*xx/sizeImageKernel).^2 ; 
                %dkern = -pi/sizeImageKernel*sin(pi/2*xx/sizeImageKernel).*cos(pi/2*xx/sizeImageKernel) ;
        end
        
    % IMAGE SMOOTHING
        KERN = kern*kern' ;
        KERN = KERN/sum(KERN(:)) ;
        Smooth = @(img)conv2(mat2gray(img),KERN,'same') ;
        
    % IMAGE DERIVATION
        if 1 % Using sparse matrices (could be extended to arbitrary kernel)
            % Indices preparation
                iii = [1:numel(indDOMAIN) 1:numel(indDOMAIN)]' ;
                jjj_dx = [sub2ind([nI nJ],IId,JJd-1) ; sub2ind([nI nJ],IId,JJd+1)] ;
                jjj_dy = [sub2ind([nI nJ],IId-1,JJd) ; sub2ind([nI nJ],IId+1,JJd)] ;
                val_deriv = 1/2*[-ones(size(indDOMAIN)) ; ones(size(indDOMAIN))] ;
            % Sparse derivation matrices
                Dx = sparse(iii,jjj_dx,val_deriv,numel(indDOMAIN),nI*nJ) ;
                Dy = sparse(iii,jjj_dy,val_deriv,numel(indDOMAIN),nI*nJ) ;
            % Restriction to the only needed domain
                gradDOMAIN = find( any(Dx,1) | any(Dy,1) )' ;
                Dx = Dx(:,gradDOMAIN) ;
                Dy = Dy(:,gradDOMAIN) ;
            % Function declaration
                dI_dx = @(img)Dx*img(gradDOMAIN) ;
                dI_dy = @(img)Dy*img(gradDOMAIN) ;
        else % Using convolution kernels
            % Apply finite differences anyway
                kern = [0 1 0]' ; 
                dkern = [1 0 -1]'/2 ;
            % 2D Kernels
                dKERN_dx = kern*dkern' ;
                dKERN_dy = dkern*kern' ;
            % Derivation function
                dI_dx = @(img)conv2(img,dKERN_dx,'same') ;
                dI_dy = @(img)conv2(img,dKERN_dy,'same') ;
        end
        
    % Transfer function to convert vector to sparse diagonal (x10 speed for big DIC meshes, welcome to Matlab...)
        spdiag = @(data) sparse(1:numel(data),1:numel(data),data(:)) ;


