%% RUN GLOBAL DIC !

lastPlotTime = tic ;
for ii = dicFrames
    % Init
        it = 0 ;
        outFlag = false ;
        dUo = ones(nNodes,2)*0 ;
    % Import and display current image
        img2i = IMG(:,:,:,ii) ;
        img2 = Smooth(img2i) ;
    % Init Valid Geometry for this frame
        VALID.Nodes(:,ii) = VALID.Nodes(:,ii-dicDir) ;
        VALID.Elems(:,ii) = VALID.Elems(:,ii-dicDir) ;
        VALID.Edges(:,ii) = VALID.Edges(:,ii-dicDir) ;
        VALID.NakedEdges(:,ii) = VALID.NakedEdges(:,ii-dicDir) ;
    % First Guess for the positions
        % Xn(:,:,ii) has already been initialized to Xn0(:,:,ii) (navDIC positions)
        % Take the previously computed frame if needed
            if ~useNavDICXn(ii) || all(isnan(Xn(:,1,ii)))
                Xn(:,:,ii) = Xn(:,:,ii-dicDir) ;
            end
        % Add the previous "correction" as convergence help
            if addPreviousCorrection && abs(refFrame-ii)>=2
                if useNavDICXn(ii) && ~all(isnan(Xn0(:,1,ii-dicDir))) % Add the correction of the previous frame with regard to the navDIC positions
                    correctionXn = Xn(:,:,ii-dicDir) - Xn0(:,:,ii-dicDir) ;
                else % Add the correction of the previous frame with regard to the before-the-previous frame
                    correctionXn = (Xn(:,:,ii-dicDir)-Xn(:,:,ii-2*dicDir)) * (frames(ii)-frames(ii-dicDir))/(frames(ii-dicDir)-frames(ii-2*dicDir)) ;
                end
                Xn(:,:,ii) = Xn(:,:,ii) + correctionXn * 1.0 ;
            end
    % Displacement guess
        Un(:,:,ii) = Xn(:,:,ii) - Nodes ;
    % Newton-Raphson
        while ~any(outFlag) && ~stopBtn.Value
            % VALID GEOMETRY 
                % OUT-OF-FRAME POINTS
                    if cullOutOfFrame
                        VALID.Nodes(:,ii) = VALID.Nodes(:,ii) & ceil(Xn(:,1,ii))<nJ+1 & floor(Xn(:,1,ii))>=1 & ceil(Xn(:,2,ii))<nI+1 & floor(Xn(:,2,ii))>=1 ;
                    end
                % Elements
                    VALID.Elems(:,ii) = VALID.Elems(:,ii) & sum(tri2nod(VALID.Nodes(:,ii),:),1)'==3 ; % triangles with still their three nodes valid
                % Edges
                    VALID.Edges(:,ii) = VALID.Edges(:,ii) & sum(edg2nod(VALID.Nodes(:,ii),:),1)'==2 ; % edges with still their two ends points valid
                    VALID.NakedEdges(:,ii) = sum(tri2edg(:,VALID.Elems(:,ii)),2)<2 ; % edges with less than two elements valid
                % Nodes
                    VALID.Nodes(:,ii) = VALID.Nodes(:,ii) & sum(tri2nod(:,VALID.Elems(:,ii)),2)>0 ; % nodes must be linked to at least one valid element
                    VALID.Nodes(:,ii) = VALID.Nodes(:,ii) & sum(edg2nod(:,VALID.Edges(:,ii)),2)>0 ; % nodes must be linked to at least one valid edge
                    nVALID = sum(VALID.Nodes(:,ii)) ;
                    if nVALID==0 ; break ; end
                    tri2nod_it = tri2nod(VALID.Nodes(:,ii),VALID.Elems(:,ii)) ;
                    valance_it = tri2nod_it./sum(tri2nod_it,2) ;
                % Valid domain at the pixel level
                    validDomain_it = any(INSIDE(:,VALID.Elems(:,ii)),2) ;
            % PIXEL-WISE DISPLACEMENT AND DIC DOMAIN
                % VALID DIC Domain
                % Compute the displacement 
                    Ui = Un(:,:,ii) ; 
                    Ui(~VALID.Nodes(:,ii),:) = 0 ;
                    Up = MAPPING*Ui ;
                % New position of each pixel
                    JJp = JJd+Up(:,1) ;
                    IIp = IId+Up(:,2) ;
                % Cull out-of-frame pixels
                    if ~cullOutOfFrame % pixels can be outside only if Nodes are outside
%                         outOfFrame = IIp<1 | IIp>nI | JJp<1 | JJp>nJ ;
                    end
            % IMAGE WARPING
                % Reduce the interpolation frame
                    iii = floor(min(IIp)):ceil(max(IIp)) ;
                    jjj = floor(min(JJp)):ceil(max(JJp)) ;
                % Warp the image
                    img2v = interp2(JJ(iii,jjj),II(iii,jjj),img2(iii,jjj),JJp,IIp,imWarpInterpOrder) ;
            % IMAGE MOMENTS
                % Select the part of the reference image of interest
                    img1v = img1(indDOMAIN) ;
                % Mean over elements
                    if refImageChanged ; meanImg1 = (WEIGHT'*img1v)./sumWEIGHT ; end
                    meanImg2 = (WEIGHT'*img2v)./sumWEIGHT ;
                % Zero-local-mean on pixels
                    if refImageChanged ; img1m = img1v-WEIGHT*meanImg1 ; end
                    img2m = img2v-WEIGHT*meanImg2 ;
                % Norm over element
                    if refImageChanged ; normImg1 = sqrt(WEIGHT'*(img1m.^2)) ; end
                    normImg2 = sqrt(WEIGHT'*(img2m.^2)) ;
                % Zero-local-mean-normalized images
                    if refImageChanged ; img1mz = img1m./(WEIGHT*normImg1) ; end
                    img2mz = img2m./(WEIGHT*normImg2) ;
            % IMAGE FUNCTIONAL
                switch diffCriterion
                    case 'Diff' % Simple difference
                        diffImg = img1v-img2v ;
                        weight = speye(2*nVALID) ;
                    case 'ZM_Diff' % Zero-mean difference
                        diffImg = img1m-img2m ;
                        weight = speye(2*nVALID) ;
                    case 'ZM_N_Diff' % Normalized Zero-mean difference
                        diffImg = img1mz-img2mz ;
                        switch size(WEIGHT,2)
                            case nNodes
                                ww = 1./normImg1(VALID.Nodes(:,ii)) ;
                            case nElems
                                ww = valance_it*normImg1(VALID.Elems(:,ii)) ;
                        end
                        weight = sparse(1:2*nVALID,1:2*nVALID,[ww;ww]) ;
                end
            % DESCENT COMPUTATION
                % Recompute Jacobian and Hessian if needed
                    if refImageChanged || strcmp(method,'full-GN')
                        % Method-dependent image used for image gradients
                            switch method
                                case 'full-GN'
                                    if it==0 % the image2 has changed
                                        dImg2_dx = conv2(img2,[1 0 -1]/2,'same') ;
                                        dImg2_dy = conv2(img2,[1 0 -1]'/2,'same') ;
                                    end
                                    dImg_dx = interp2(JJ(iii,jjj),II(iii,jjj),dImg2_dx(iii,jjj),JJp,IIp,imWarpInterpOrder) ;
                                    dImg_dy = interp2(JJ(iii,jjj),II(iii,jjj),dImg2_dy(iii,jjj),JJp,IIp,imWarpInterpOrder) ;
                                otherwise
                                    dImg_dx = dI_dx(img1) ;
                                    dImg_dy = dI_dy(img1) ;
                            end
                        % Jacobian
                            dImg_da = [...
                                    spdiag(dImg_dx)*MAPPING ... 
                                    spdiag(dImg_dy)*MAPPING ... 
                                    ] ;
                        % Hessian
                            Hess = dImg_da'*dImg_da ;
                        % No need to compute next time
                            refImageChanged = false ;
                    end
                % Compute the first RMSE derivative
                    dr_da = (dImg_da(validDomain_it,:)'*diffImg(validDomain_it)) ;
                % Contraint on the SECOND displacement gradient
                    validDOF = [VALID.Nodes(:,ii);VALID.Nodes(:,ii)] ;
                    wCon = ones(1,nVALID) ;
                    switch strainCriterion
                        case 'normal'
                            vEdg = repmat(~VALID.NakedEdges(:,ii),[2 1]) ;
                            vEle = [VALID.Elems(:,ii);VALID.Elems(:,ii);VALID.Elems(:,ii)] ;
                            CONS = B(vEle,validDOF)'*(Ed(vEdg,vEle)'*Ed(vEdg,vEle))*B(vEle,validDOF) ;
                        case 'full'
                            vEdg = repmat(~VALID.NakedEdges(:,ii),[4 1]) ;
                            vEle = [VALID.Elems(:,ii);VALID.Elems(:,ii);VALID.Elems(:,ii);VALID.Elems(:,ii)] ;
                            CONS = G(vEle,validDOF)'*(Ed(vEdg,vEle)'*Ed(vEdg,vEle))*G(vEle,validDOF) ;
                            %if strcmp(regCrit,'rel') ; wCon = 1./max(epsTrsh,abs(repmat(valance_it,[1 4])*G(vEle,validDOF)*[Un(VALID.Nodes(:,ii),1,ii);Un(VALID.Nodes(:,ii),2,ii)])) ; end
                    end
                    wCon = sparse(1:2*nVALID,1:2*nVALID,[wCon;wCon]) ;
                    DU = reshape(Un(VALID.Nodes(:,ii),:,ii),[],1) ;
                    switch regCrit
                        case 'abs' % Strain computed from the reference state
                        case 'rel' % Strain computed from the previous state
                            if abs(refFrame-ii)>=2
                                DU = DU - reshape(Un(VALID.Nodes(:,ii),:,ii-dicDir),[],1) ;
                            end
                    end
                % Updating DOFs, X*a=b
                    X = ( ...
                            weight*Hess(validDOF,validDOF)*weight...
                            + beta*wCon*CONS*wCon...
                        ) ; 
                    b = (...
                            weight*dr_da(validDOF)...
                            - beta*wCon*CONS*DU...
                        ) ;
                    a = X\b ;
                    descent = -b'*a ; % Descent direction (should be always<0)
                % Displacement
                    dU = reshape(a,[nVALID 2]) ;
                % Descent Step
                    % Displacement oscillations
                        corr_dU = sum(dUo(validDOF).*dU(:))/(norm(dUo(validDOF))*norm(dU(:))) ;
                        dUo(:) = 0 ; dUo(validDOF) = dU(:) ;
                    % Compute the step
                        step = stepRatio ;
                        if corr_dU>.95
                            %step = 1-(1-stepRatio)^4 ;
                        end
                        if corr_dU<-.1 % Need to decrease the step
                            step = stepRatio^(1/(1+corr_dU)) ;
                        end
                % Positions
                    Un(VALID.Nodes(:,ii),:,ii) = Un(VALID.Nodes(:,ii),:,ii) + dU * step ;
                    Xn(:,:,ii) = Nodes + Un(:,:,ii) ;
            % CONVERGENCE CITERIONS
                % number of iterations
                    it = it+1 ;
                % Residues Maps
                    RMSE(ii,it) = norm(X*a-b) ; %sqrt(sum(diffImg(:).^2)) ;
                    relativeResidueVariation = (RMSE(ii,it)-RMSE(ii,max(1,it-1)))/RMSE(ii,max(1,it-1)) ;
                    meanSquaredElemResid = ((WEIGHT'*abs(diffImg))) ;
                    corrCoeff = abs(WEIGHT'*(img1mz.*img2mz)) ;
                % Criterions
                    normA = max(abs(a)) ; norm(a)/nNodes ;
                % Convergence criterion
                    if it>maxIt ; outFlag = 'number of iterations' ; end
                    if normA<minNorm ; outFlag = 'norm of descent direction' ; end
                    if corr_dU<minCorrdU ; outFlag = 'descent direction oscillates' ; end
                    if it>1 && relativeResidueVariation>maxResidueRelativeVariation ; outFlag = 'RMSE did not decreased sufficiently' ; end
            % POINTS/ELEMENTS VALIDATION/DELETION
                % DECORRELATED ELEMENTS
                    cullGeo = corrCoeff<minCorrCoeff | meanSquaredElemResid>maxMeanElemResidue ;
                    if any(cullGeo) && (any(outFlag) || (normA/minNorm)<thresholdValidGeometry)
                        switch size(WEIGHT,2) 
                            case nElems % Correlation at the element level
                                VALID.Elems(VALID.Elems(:,ii),ii) = VALID.Elems(VALID.Elems(:,ii),ii) & ~cullGeo ; 
                            case nNodes % Correlation at the node level$
                                VALID.Nodes(VALID.Nodes(:,ii),ii) = VALID.Nodes(VALID.Nodes(:,ii),ii) & cullGeo ;
                        end
                        outFlag = false ; % Force an new iteration
                    end
                % SET THE NON-VALID NODES TO NAN
                    Xn(~VALID.Nodes(:,ii),:,ii) = NaN ;
                    Un(~VALID.Nodes(:,ii),:,ii) = NaN ;
            % DISPLAY
                % Always display infos
                    infosText.String = [ ' INFOS' ...
                                         ' | Frame: ' num2str(frames(ii)) ...
                                         ' | It: ' num2str(it) ...
                                         ' | norm.A: ' num2str(normA,2) ...
                                         ' | corr.dU: ' num2str(corr_dU,3) ...
                                         ' | RMSE: ' num2str(RMSE(ii,it),5) ...
                                         ' (' num2str(relativeResidueVariation,3) ')' ...
                                         ' | outFlag: ' mat2str(outFlag) ...
                                         ] ; 
                % Heavy Plots
                    if plotEachIteration || (any(outFlag) && plotEachFrame) || toc(lastPlotTime)>1/plotRate % (any(outFlag) && toc(lastPlotTime)>1/plotRate)
                        im.CData = repmat(img2,[1 1 3]) ;
                        %%
                        residues = ... WEIGHT*corrCoeff ... Correlation coeffient
                                    abs(diffImg) ... NR residues
                                   ... Up(:,1) ... Ux displacement
                                   ... WEIGHT*meanSquaredElemResid ... Sum of squared NR residues
                                   ... WEIGHT*((WEIGHT'*abs(diffImg))./sumWEIGHT) ...
                                   ;
                        imRes.CData(indDOMAIN) = residues ;
                        %%
                        mesh.Vertices = Xn(:,:,ii) ;
                        mesh.Faces = Elems(VALID.Elems(:,ii),:) ;
                        mesh.FaceVertexCData = ones(nNodes,1)*NaN ; mesh.FaceVertexCData(VALID.Nodes(:,ii)) = sqrt(sum(dU.^2,2)) ;
                        %figure(figDebug) ; clf ; ind = (0:nJ-1)*nI+ceil(nI/2) ; plot(img1v(ind)) ; plot(img2v(ind)) ;
                        lastPlotTime = tic ; 
                    end
            % Pause execution ?
                if pauseAtPlot && ~continueBtn.Value
                    while ~stopBtn.Value && ~nextBtn.Value && ~continueBtn.Value
                        drawnow ;
                    end
                    nextBtn.Value = 0 ;
                    if continueBtn.Value
                        nextBtn.Visible = 'off' ;
                        continueBtn.Visible = 'off' ;
                    end
                end
            % Draw
                drawnow ;
        end
        % Modify the reference image if needed
            if weightCurrentImage>0
                img1(indDOMAIN) = img1(indDOMAIN)*(1-weightCurrentImage) + img2v*weightCurrentImage ;
                refImageChanged = true ;
            end
        % Out Criterions
            if stopBtn.Value; break ; end
            if nVALID==0; break ; end
end