function MovingPoints = fftDispMethod(PtsMov,PtsRef,imgMov,imgRef,CorrSize)

    % PARAMETERS
        dir = 'both' ; % displacement directions: 'both', 'X' or 'Y'
        CorrSize = 81*[1 1] ; % Rectangular window
        m = round(CorrSize/4) ; % Margin to truncate borders
        uMax = CorrSize/2 ; 30*[1 1] ; % Maximum allowed displacement per iteration
        iterateIfMaxDispHigherThan = 1 ;
        maxIt = 10 ;
        FIT =   ... 'LS' ...
                 'SVD' ...
                ; 
        windowing = true ; % apply a blackman windowing
        
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
                        [JJ,II] = meshgrid(jj,ii) ;
                        III = repmat(II,[1 1 nPts]) ;
                        JJJ = repmat(JJ,[1 1 nPts]) ;
                    % Imagette for each point
                        if it==1 ; IIIRef = bsxfun(@plus,III,reshape(disPtsRef(:,2),[1 1 nPts])) ; end
                        if it==1 ; JJJRef = bsxfun(@plus,JJJ,reshape(disPtsRef(:,1),[1 1 nPts])) ; end
                        IIIMov = bsxfun(@plus,III,reshape(disPtsMov(:,2),[1 1 nPts])) ;
                        JJJMov = bsxfun(@plus,JJJ,reshape(disPtsMov(:,1),[1 1 nPts])) ;
                    % If outside of Image
                        Valid = IIIRef>0 ...
                                & IIIRef<=imgSize(1) ...
                                & JJJRef>0 ...
                                & JJJRef<=imgSize(2) ...
                                & IIIMov>0 ...
                                & IIIMov<=imgSize(1) ...
                                & JJJMov>0 ...
                                & JJJMov<=imgSize(2) ;
                        Valid = ~any(any(~Valid,1),2) ;
                        %notValid = ~Valid ;
                        nPtsValid = sum(Valid(:)) ;
                    % Linear Indices
                        NNNRef = sub2ind(size(imgRef),IIIRef(:,:,Valid),JJJRef(:,:,Valid)) ;
                        NNNMov = sub2ind(size(imgMov),IIIMov(:,:,Valid),JJJMov(:,:,Valid)) ;
                    % Pixel Values
                        imagettesRef = imgRef(NNNRef) ;
                        imagettesMov = imgMov(NNNMov) ;

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
                    N = floor(CorrSize/2) ;
                    %H = hankel(0:N(1)-1,N(1)-1:CorrSize(1)-1),CorrSize(1)*hankel(
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
                            % SVD Filtering
                                if 1
                                    [U,w] = eig(phi*phi','vector') ;
                                    [w,ind] = sort(w,'descend') ;
                                    U = U(:,ind(1)) ;
                                    V = phi'*diag(1./sqrt(w(1)))*U ;
                                else
                                    %[U,~,V] = svd(phi,'econ') ;
                                end
                            % Shift Invariance
                                k(p,2) = U(1:end-1,1)\U(2:end,1) ;
                                k(p,1) = conj(V(1:end-1,1)\V(2:end,1)) ;
                        end
                end


            % Displacement
                U = real(1i*log(k)) ;
                U = bsxfun(@times,U,flip(CorrSize(:)')/2/pi);
                normU = sqrt(sum(U.^2,2)) ;
                U(normU>uMax(2),1) = NaN ;
                U(normU>uMax(1),2) = NaN ;

            % Zero components if needed
                switch dir
                    case 'both'
                    case 'X'
                        U(:,2) = 0 ; % Uy = 0
                    case 'Y'
                        U(:,1) = 0 ; % Ux = 0
                end

            % Moving Points
                MovingPoints(Valid,:) = U + disPtsMov(Valid,:) - (disPtsRef(Valid,:)-PtsRef(Valid,:)) ;

            % Absolute value of the displacement increment
                normU = zeros(nPts,1) ;
                normU(Valid) = sqrt(sum((MovingPoints(Valid,:)-PtsMov(Valid,:)).^2,2)) ;
                
            % Apply the increment to PtsMov
                PtsMov(Valid,:) = MovingPoints(Valid,:) ;
                
            % Delete PtsMov for which the increment was smaller than the
            % criterion
                PtsMov(normU<iterateIfMaxDispHigherThan,:) = NaN ;
                if all(isnan(PtsMov(:))) ; outFlag = true ; end
                
            
            if it>=maxIt ; outFlag = true ; end   
            
        end
        
    % Timing
        disp(['fftdisp: ' num2str(toc(startTime)*1000,'%.1f') ' ms (' num2str(it) ' iteration(s))']) ;  
        
end

        