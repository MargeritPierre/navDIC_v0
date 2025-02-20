function MovingPoints = fftDispMethod(PtsMov,PtsRef,imgMov,imgRef,params)

    % PARAMETERS
        dir = 'both' ; % displacement directions: 'both', 'X' or 'Y'
        CorrSize = 301*[1 1] ; % Rectangular window
        m = round(CorrSize/4) ; % Margin to truncate borders
        uMax = CorrSize/2 ; 30*[1 1] ; % Maximum allowed displacement per iteration
        maxImagetShift = CorrSize/2 ; % Maximum allowed imagette shift (close to image borders)
        iterateIfMaxDispHigherThan = 1 ;
        maxIt = 10 ;
        FIT =   ...'LS' ...
                ... 'SVD' ...
                 'COR' ...
                ; 
        spatialSmoothing = prod(CorrSize)<20000 ;
        spatialSmoothingRatio = 1/2 ;
        windowing = true ; % apply a blackman windowing
    
    % Override parameters (not very nice..)
        if ~isempty(params)
            for ff = fieldnames(params)'
                eval([ff{1} ' = params.' ff{1} ';' ]) ;
            end
        end
        
    % INITIALIZE
        MovingPoints = zeros(size(PtsMov))*NaN ;
        
    % MODIY PARAMETERS
        switch dir
            case 'X'
                m(1) = 0 ;
            case 'Y'
                m(2) = 0 ;
        end
    
    % INFOS
        nPts = size(PtsMov,1) ;
        imgSize = size(imgRef) ;
        startTime = tic ;
        
        
    % LOOP !!!
        it = 0 ;
        outFlag = false ;
        while ~outFlag
            
            it = it+1 ;
            
            % Discrete positions
                if it==1 ; disPtsRef = round(PtsRef) ; end
                disPtsMov = round(PtsMov) ;

            % Create Imagettes
                % Indices
                    % One imagette centered in 0
                        ii = floor(-CorrSize(1)/2+(1:CorrSize(1))) ; % -floor(CorrSize(1)/2)+(1:CorrSize(1)) ;
                        jj = floor(-CorrSize(2)/2+(1:CorrSize(2))) ; % -floor(CorrSize(2)/2)+(1:CorrSize(2)) ;
                    % All imagettes centerd in 0
                        [JJ,II] = meshgrid(jj,ii) ;
                        III = repmat(II,[1 1 nPts]) ;
                        JJJ = repmat(JJ,[1 1 nPts]) ;
                    % Imagette for each point
%                         if it==1 ; IIIRef = bsxfun(@plus,III,reshape(disPtsRef(:,2),[1 1 nPts])) ; end
%                         if it==1 ; JJJRef = bsxfun(@plus,JJJ,reshape(disPtsRef(:,1),[1 1 nPts])) ; end
                        IIIRef = bsxfun(@plus,III,reshape(disPtsRef(:,2),[1 1 nPts])) ;
                        JJJRef = bsxfun(@plus,JJJ,reshape(disPtsRef(:,1),[1 1 nPts])) ;
                        IIIMov = bsxfun(@plus,III,reshape(disPtsMov(:,2),[1 1 nPts])) ;
                        JJJMov = bsxfun(@plus,JJJ,reshape(disPtsMov(:,1),[1 1 nPts])) ;
                    % Shift imagettes if too close from boundaries
                        % Get the overlapping pixels
                            onTop = min(min(IIIRef(1,1,:)-1,0),min(IIIMov(1,1,:)-1,0)) ;
                            onBottom = max(max(IIIRef(end,end,:)-imgSize(1),0),max(IIIMov(end,end,:)-imgSize(1),0)) ;
                            onLeft = min(min(JJJRef(1,1,:)-1,0),min(JJJMov(1,1,:)-1,0)) ;
                            onRight = max(max(JJJRef(end,end,:)-imgSize(2),0),max(JJJMov(end,end,:)-imgSize(2),0)) ;
                        % Compute the shifts
                            imagetShiftII = onTop ; 
                            imagetShiftII(logical(onBottom)) = onBottom(logical(onBottom)) ;
                            imagetShiftII(logical(onBottom & onTop)) = NaN ;
                            imagetShiftJJ = onLeft ; 
                            imagetShiftJJ(logical(onRight)) = onRight(logical(onRight)) ;
                            imagetShiftJJ(logical(onLeft & onRight)) = NaN ;
                        % Prevent a too big shift
                            imagetShiftII(abs(imagetShiftII)>maxImagetShift(1)) = NaN ;
                            imagetShiftJJ(abs(imagetShiftJJ)>maxImagetShift(2)) = NaN ;
                        % Apply
                            IIIMov = IIIMov - imagetShiftII ;
                            JJJMov = JJJMov - imagetShiftJJ ;
                            IIIRef = IIIRef - imagetShiftII ;
                            JJJRef = JJJRef - imagetShiftJJ ;
                    % If points outside of Image
                        ValidPts = ~isnan(IIIRef(1,1,:) + JJJRef(1,1,:) + IIIMov(1,1,:) + JJJMov(1,1,:)) ;
%                         ValidIJ = IIIRef>0 ...
%                                 & IIIRef<=imgSize(1) ...
%                                 & JJJRef>0 ...
%                                 & JJJRef<=imgSize(2) ...
%                                 & IIIMov>0 ...
%                                 & IIIMov<=imgSize(1) ...
%                                 & JJJMov>0 ...
%                                 & JJJMov<=imgSize(2) ;
%                         ValidPts = ~any(any(~ValidIJ,1),2) ;
                        %notValid = ~Valid ;
                        nPtsValid = sum(ValidPts(:)) ;
                    % Linear Indices
                        NNNRef = sub2ind(size(imgRef),IIIRef(:,:,ValidPts),JJJRef(:,:,ValidPts)) ;
                        NNNMov = sub2ind(size(imgMov),IIIMov(:,:,ValidPts),JJJMov(:,:,ValidPts)) ;
                    % Pixel Values
                        imagettesRef = imgRef(NNNRef) ;
                        imagettesMov = imgMov(NNNMov) ;
                        
            % Switch to double data
                imagettesRef = double(imagettesRef) ;
                imagettesMov = double(imagettesMov) ;
                        
            % Zero-mean
                imagettesRef = imagettesRef-reshape(mean(reshape(imagettesRef,[],nPtsValid),1),[1 1 nPtsValid]) ;
                imagettesMov = imagettesMov-reshape(mean(reshape(imagettesMov,[],nPtsValid),1),[1 1 nPtsValid]) ;

            % Windowing
                if windowing
                    switch dir
                        case 'both'
                            WIN = blackman(CorrSize(1))*blackman(CorrSize(2))' ;
                        case 'X'
                            WIN = blackman(CorrSize(1))*ones(1,CorrSize(2)) ;
                        case 'Y'
                            WIN = ones(CorrSize(1),1)*blackman(CorrSize(2))' ;
                    end
                    WIN = WIN./sum(WIN(:)) ;
                    imagettesRef = bsxfun(@times,double(imagettesRef),WIN) ;
                    imagettesMov = bsxfun(@times,double(imagettesMov),WIN) ;
                end

            % FFT
                switch dir
                    case 'both'
                        fftRef = fft(fft(imagettesRef,[],1),[],2) ;
                        fftMov = fft(fft(imagettesMov,[],1),[],2) ;
                    case 'X'
                        fftRef = fft(imagettesRef,[],2) ;
                        fftMov = fft(imagettesMov,[],2) ;
                    case 'Y'
                        fftRef = fft(imagettesRef,[],1) ;
                        fftMov = fft(imagettesMov,[],1) ;
                end

            % Phase field
                PHI = fftMov.*conj(fftRef) ;
                PHI = fftshift(fftshift(PHI,1),2) ;
                PHI = PHI./abs(PHI+eps) ;

            % Hankel Matrix
                if strcmp(FIT,'SVD')
                    indI = 1+m(1):CorrSize(1)-m(1) ; 
                    indJ = 1+m(2):CorrSize(2)-m(2) ;
                    if spatialSmoothing
                        LL = [numel(indI) numel(indJ)] ;
                        KK = floor(LL*spatialSmoothingRatio) ;
                        H1 = 1 + repmat( hankel(0:KK(1)-1,KK(1)-1:LL(1)-1) , [KK(2) LL(2)-KK(2)+1]) ;
                        H2 = 1 + kron(hankel(0:KK(2)-1,KK(2)-1:LL(2)-1) ,ones(KK(1),LL(1)-KK(1)+1)) ;
                        HH = H1 + LL(1)*(H2-1) ;
                    end
                end


            % Wavevectors
                k = zeros(nPtsValid,2) ;
                switch FIT
                    case 'LS' % Fast implementation
                        % Reshape the data
                            PHI = reshape(PHI,[],nPtsValid) ;
                        % Selection indices
                            indUpI = [false(m(1),CorrSize(2)) ; ...
                                        [false(CorrSize(1)-2*m(1)-1,m(2)),...
                                            true(CorrSize(1)-2*m(1)-1,CorrSize(2)-2*m(2)),...
                                            false(CorrSize(1)-2*m(1)-1,m(2))] ; ...
                                       false(m(1)+1,CorrSize(2))] ; 
                            indDwnI = [false(m(1)+1,CorrSize(2)) ; ...
                                        [false(CorrSize(1)-2*m(1)-1,m(2)),...
                                            true(CorrSize(1)-2*m(1)-1,CorrSize(2)-2*m(2)),...
                                            false(CorrSize(1)-2*m(1)-1,m(2))] ; ...
                                       false(m(1),CorrSize(2))] ; 
                            indUpJ = [false(CorrSize(1),m(2)) , ...
                                        [false(m(1),CorrSize(2)-2*m(2)-1); ...
                                            true(CorrSize(1)-2*m(1),CorrSize(2)-2*m(2)-1);...
                                            false(m(1),CorrSize(2)-2*m(2)-1)] , ...
                                       false(CorrSize(1),m(2)+1)] ; 
                            indDwnJ = [false(CorrSize(1),m(2)+1) , ...
                                        [false(m(1),CorrSize(2)-2*m(2)-1); ...
                                            true(CorrSize(1)-2*m(1),CorrSize(2)-2*m(2)-1);...
                                            false(m(1),CorrSize(2)-2*m(2)-1)] , ...
                                       false(CorrSize(1),m(2))] ; 
                        % Shifted data
                            phiUpI = PHI(indUpI(:),:) ;
                            phiDwnI = PHI(indDwnI(:),:) ;
                            phiUpJ = PHI(indUpJ(:),:) ;
                            phiDwnJ = PHI(indDwnJ(:),:) ;
                        % LS Formatting: u1\u2 = (u1'*u2)/|u1|^2
                            normPhiUpI = sum(abs(phiUpI).^2,1) ;
                            normPhiUpJ = sum(abs(phiUpJ).^2,1) ;
                            scalPhiI = sum(conj(phiUpI).*phiDwnI,1) ;
                            scalPhiJ = sum(conj(phiUpJ).*phiDwnJ,1) ;
                        % LS Estimation
                            k(:,2) = scalPhiI./normPhiUpI ;
                            k(:,1) = scalPhiJ./normPhiUpJ ;
                    case 'SVD' % Use SVD Decomposition (Slower)
                        for p = 1:nPtsValid
                            % Select and truncate
                                phi = PHI(indI,indJ,p) ;
                            % Spatial smoohing
                                if spatialSmoothing ; phi = phi(HH) ; end
                            % SVD Filtering
                                Cpp = phi*phi' ;
                                [W,w] = eig(Cpp,'vector') ;
                                [~,ind] = sort(w,'descend') ;
                                W = W(:,ind(1)) ;
                            % Shift Invariance
                                if spatialSmoothing
                                    k(p,2) = W(H1(:,1)<KK(1),1)\W(H1(:,1)>1,1) ;
                                    k(p,1) = W(H2(:,1)<KK(2),1)\W(H2(:,1)>1,1) ;
                                else
                                    W = phi'*diag(1./sqrt(w(1)))*W ;
                                    k(p,2) = W(1:end-1,1)\W(2:end,1) ;
                                    k(p,1) = conj(W(1:end-1,1)\W(2:end,1)) ;
                                end
                        end
                    case 'COR'
                    % kept indices
                        indI = 1+m(1):CorrSize(1)-m(1) ; 
                        indJ = 1+m(2):CorrSize(2)-m(2) ;
                        phi = PHI(indI,indJ,:) ;
                    % Sizes
                        L = size(phi,[1 2]) ; N = ceil(L/2) ; M = L-N+1 ;
                    % Fourier transform pre-compute
                        Fphi = fft(fft(phi,[],1),[],2) ;
                    % Intialize
                        W = ones([M nPtsValid])/sqrt(prod(M)) ;
                        dW = inf ;
                    % Eigenvector estimation via Power iterations w = (Cpp*w)/norm(Cpp*w)
                        while max(abs(dW(:)))>1e-6
                        % product of Cpp = Hp*Hp' with W
                            W0 = W ;
                        % W = Hp'*W
                            W = flip(flip(conj(W),1),2) ;
                            Fw = fft(fft(W,L(1),1),L(2),2) ;
                            Fw = Fw.*Fphi ;
                            W = ifft(ifft(Fw,[],1),[],2) ;
                            W = W(M(1):L(1),M(2):L(2),:) ;
                            W = conj(W) ;
                        % W = Hp*W ;
                            W = flip(flip(W,1),2) ;
                            Fw = fft(fft(W,L(1),1),L(2),2) ;
                            Fw = Fw.*Fphi ;
                            W = ifft(ifft(Fw,[],1),[],2) ;
                            W = W(N(1):L(1),N(2):L(2),:) ;
                        % W = W/norm(W)
                            W = W./sqrt(sum(abs(W).^2,[1 2])) ;
                        % Update
                            dW = W-W0 ;
                        end
                    % Least-squares pole estimation
                        k(:,2) = sum(conj(W(1:end-1,:,:)).*W(2:end,:,:),[1 2])./sum(abs(W(1:end-1,:,:)).^2,[1 2]) ;
                        k(:,1) = sum(conj(W(:,1:end-1,:)).*W(:,2:end,:),[1 2])./sum(abs(W(:,1:end-1,:)).^2,[1 2]) ;
                end


            % Displacement
                W = real(1i*log(k)) ;
                W = bsxfun(@times,W,flip(CorrSize(:)')/2/pi);
                normU = sqrt(sum(W.^2,2)) ;
                W(normU>uMax(2),1) = NaN ;
                W(normU>uMax(1),2) = NaN ;

            % Zero components if needed
                switch dir
                    case 'both'
                    case 'X'
                        W(:,2) = 0 ; % Uy = 0
                    case 'Y'
                        W(:,1) = 0 ; % Ux = 0
                end

            % Moving Points
                MovingPoints(ValidPts,:) = W + disPtsMov(ValidPts,:) - (disPtsRef(ValidPts,:)-PtsRef(ValidPts,:)) ;

            % Absolute value of the displacement increment
                normU = zeros(nPts,1) ;
                normU(ValidPts) = sqrt(sum((MovingPoints(ValidPts,:)-PtsMov(ValidPts,:)).^2,2)) ;
                
            % Apply the increment to PtsMov
                PtsMov(ValidPts,:) = MovingPoints(ValidPts,:) ;
                
            % Delete PtsMov for which the increment was smaller than the
            % criterion
                PtsMov(normU<iterateIfMaxDispHigherThan,:) = NaN ;
                if all(isnan(PtsMov(:))) ; outFlag = true ; end
                
            
            if it>=maxIt ; outFlag = true ; end   
            
        end
        
    % Timing
        disp(['fftdisp: ' num2str(toc(startTime)*1000,'%.1f') ' ms (' num2str(it) ' iteration(s))']) ;  
        
end

        