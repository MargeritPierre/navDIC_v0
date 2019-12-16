function MovingPoints = icgnCorrMethod(PtsMov,PtsRef,imgMov,imgRef,CorrSize)
    
    % Convertir image en 8 bits en float 
    t = tic ;
    tic
        if isa(imgMov(1,1),'uint8')
            imgMov = double(imgMov)/255 ;
        end
        if isa(imgRef(1,1),'uint8')
            imgRef = double(imgRef)/255 ;
        end
    time = toc;
    disp(['fin conversion : ', num2str(time), ' s.']) ; 
    tic ;
    % PARAMETERS
        dir = 'both' ; % displacement directions: 'both', 'X' or 'Y'
        transf = 'cauchy' ; %'uniaxiale' ; % false ; %
        %local = true ; % 'transformation sur une imagette (sinon transformation de l'image complete)'
        interpolMethod = 'linear' ; %'pchip' ; %
        % carr? des Residus maximal
        R2max = 0.05 ;
        
        % fenetre de correlation
        if mod(CorrSize(1),2) == 0 
            Xfen = CorrSize(1) ;
        else
            Xfen = CorrSize(1)+1 ;
        end
        if mod(CorrSize(2),2) == 0 
            Yfen = CorrSize(2) ;
        else
            Yfen = CorrSize(2)+1 ;
        end
    
        % Definition des fonction de transformation de l'image actuelle
        % (fonction de la transformation autoris?e)

        switch transf
            case false
                uDicG = @(ptsFen,h) ptsFen + repmat(reshape(h,[1,1,2]), [Yfen+1, Xfen+1]) ;

            case 'uniaxiale'
                uDicG = @(ptsFen,ptCorr,h) repmat(reshape(h(1:2),[1,1,2]), [Yfen+1, Xfen+1]) .*...
                    (ptsFen-repmat(reshape(ptCorr,[1,1,2]),[Yfen+1, Xfen+1])) +...
                    repmat(reshape(h(3:4),[1,1,2]), [Yfen+1, Xfen+1]) +...
                    repmat(reshape(ptCorr,[1,1,2]),[Yfen+1, Xfen+1]) ;

            case 'cauchy'
                uDicG = @(ptsFen,ptCorr,h) repmat(reshape(h(1:2),[1,1,2]), [Yfen+1, Xfen+1]) .*...
                    (ptsFen-repmat(reshape(ptCorr,[1,1,2]),[Yfen+1, Xfen+1])) +...
                    repmat(reshape(h(3:4),[1,1,2]), [Yfen+1, Xfen+1]) .*...
                    (cat(3,ptsFen(:,:,2),ptsFen(:,:,1))-cat(3,repmat(ptCorr(2),[Yfen+1, Xfen+1]),repmat(ptCorr(1),[Yfen+1, Xfen+1]))) +...
                    repmat(reshape(h(5:6),[1,1,2]), [Yfen+1, Xfen+1]) +...
                    repmat(reshape(ptCorr,[1,1,2]),[Yfen+1, Xfen+1]) ;
        end
    
    % INFOS
        nPts = size(PtsMov,1) ;
        startTime = tic ;
        sizeImgs = (Xfen+1)*(Yfen+1) ; % taille imagette 
        
    % DISCRETE POSITION
        disPtsRef = round(PtsRef) ;
        time= toc ;
        disp(['temps d''initiation : ', num2str(time), ' s.']) ;
    % Initialisation moving point   
    
        MovingPoints = zeros(size(disPtsRef));
        
        for i = 1:nPts % pourrait etre remplace par une optimisation globale 
            % mais risque de ralentir si certains points diverges
            disp([ 'point n° ',num2str(i), '.']) ;
            tic
            % Deplacement a priori ( Valeur de l'image precedente )
            ui_m1 = PtsMov(i,:) - PtsRef(i,:) ;
            switch transf
                case false
                   depApriori = ui_m1' ;
                case 'uniaxiale'
                    depApriori = [1;1;ui_m1'] ;
                case 'cauchy'
                    depApriori = [1;1;0;0;ui_m1'] ;
            end
            
            % Imagettes de reference
            %----------------------------------------------------------
            % Coordonnees discretes
            imgsRef = imgRef(disPtsRef(i,2) - Yfen/2:disPtsRef(i,2) + Yfen/2,...
                disPtsRef(i,1) - Xfen/2:disPtsRef(i,1) + Xfen/2) ;
            
            % Calcul des derivees
            dFx = 0.5 * ( imgRef(disPtsRef(i,2) - Yfen/2:disPtsRef(i,2) + Yfen/2,...
                disPtsRef(i,1) - Xfen/2+1:disPtsRef(i,1) + Xfen/2+1) - ...
                imgRef(disPtsRef(i,2) - Yfen/2:disPtsRef(i,2) + Yfen/2,...
                disPtsRef(i,1) - Xfen/2-1:disPtsRef(i,1) + Xfen/2-1) ) ; 
            dFy = 0.5 * ( imgRef(disPtsRef(i,2) - Yfen/2+1:disPtsRef(i,2) + Yfen/2+1,...
                disPtsRef(i,1) - Xfen/2+1:disPtsRef(i,1) + Xfen/2+1) - ...
                imgRef(disPtsRef(i,2) - Yfen/2-1:disPtsRef(i,2) + Yfen/2-1,...
                disPtsRef(i,1) - Xfen/2-1:disPtsRef(i,1) + Xfen/2-1) ) ;
            
            % Definition de la jacobienne (fonction de la transformation autorisee)
            
            switch transf
                case false
                    % Minimisation de E = sum( f(x) - g(x+b) )
                    X = [reshape(dFx,[],1), reshape(dFy,[],1)] ;

                case 'uniaxiale'
                    % Minimisation de E = sum( f(x) - g(ax+b) ) = sum( f(x) - (g(ax+b) - aij xj dg,i - ui dg,i) 
                    % avec h = (a11; a22; u1; u2)
                    absF = repmat(-Xfen/2:Xfen/2, [Yfen+1 1]) ;
                    ordF = repmat((-Yfen/2:Yfen/2)', [1  Xfen+1 ]) ;
                    
                    X = [reshape(absF .* dFx,[],1), reshape(ordF.* dFy,[],1),...
                      reshape(dFx,[],1), reshape(dFy,[],1)] ;
                    
                case 'cauchy'
                    % Minimisation de E = sum( f(x) - g(ax+b) ) = sum( f(x) - (g(ax+b) - aij xj dg,i - ui dg,i) 
                    % avec h = (a11 ; a22 ; a12; a21; u1; u2)
                    absF = repmat(-Xfen/2:Xfen/2, [Yfen+1 1]) ;
                    ordF = repmat((-Yfen/2:Yfen/2)', [1  Xfen+1 ]) ;
                    
                    X = [reshape(absF .* dFx,[],1), reshape(ordF.* dFy,[],1), reshape(ordF.*dFx,[],1),...
                       reshape(absF .*dFy,[],1), reshape(dFx,[],1), reshape(dFy,[],1)] ;
            end
            time = toc ;
            disp( ['Calcul derivee X: ', num2str(time),' s.'])  ;
            % INITIALISATION
            h = depApriori ; % vecteur de transformation
            r = ones(length(h),1) ; % reste
            it = 0 ; % iteration
            
            % Calcul de la Hessienne approch?e (DL a l'ordre 2)
            H = (X'*X) ;
            
            % Grille imagette non transform?e configuration actuelle
            [ptsImgsRef(:,:,1),ptsImgsRef(:,:,2)] = meshgrid(disPtsRef(i,1)-Xfen/2:disPtsRef(i,1)+Xfen/2,...
               disPtsRef(i,2)-Yfen/2:disPtsRef(i,2)+Yfen/2) ;
           
            time = toc ; 
            disp( [ 'initialisation : ', num2str(time) ,' s.' ] ) ;
            tic
            while sum(r) > 10^-6 && it<50
                movPtsFen = uDicG(ptsImgsRef,disPtsRef(i,:),h) ;
                imgsMov = interp2(imgMov,movPtsFen(:,:,1),movPtsFen(:,:,2),interpolMethod) ;
                y = reshape((imgsRef - imgsMov),[],1) ;
                hcor = H \ (X'*y)   ;
                r = abs(hcor) ;
                h = h+hcor ;
                it=it+1 ;
            end
            time = toc ; 
            disp( [ 'temps d''iterations : ', num2str(time) ,' s. (' , num2str(it), ' iterations)' ] ) ;
            tic
            switch transf                
                case false
                    % Moving Points et displacement
                    if sum(r) < 10^-6
                        u = h(end-1:end)' ;
                        MovingPoints(i,:) = PtsRef(i,:) + u' ;
                    else
                        MovingPoints(i,:) = [NaN,NaN] ;
                    end
                    
                case 'uniaxiale'
                    if sum(r) < 10^-6
                        A = [h(1) 0; 0 h(2)] ;
                        %Am1 = inv(A);
                        u = -A \ h(end-1:end) + ( A \ ( PtsRef(i,:) - disPtsRef(i,:) )' - ...
                            ( PtsRef(i,:) - disPtsRef(i,:))' ) ;
                        MovingPoints(i,:) = PtsRef(i,:) - u' ;
                    else
                        MovingPoints(i,:) = [NaN,NaN] ;
                    end
                    
                case 'cauchy'
                    if sum(r) < 10^-6 && sum(sum((imgsRef - imgsMov).^2))/ sizeImgs <= R2max
                        A = [h(1) h(3); h(4) h(2)] ;
                        %Am1 = inv(A);
                        u = - A \ h(end-1:end) + ( A \ ( PtsRef(i,:) - disPtsRef(i,:) )' - ...
                            ( PtsRef(i,:) - disPtsRef(i,:))' ) ;
                        MovingPoints(i,:) = PtsRef(i,:) - u' ;
                    else
                        MovingPoints(i,:) = [NaN,NaN] ;
                    end
            end
            time = toc ; 
            disp( [ 'Temps de mise à jour resultat : ', num2str(time) ,' s.' ] ) ;
            %disp(['Coefficient de Correlation (moyenne des carres des residus) : ' ,...
            % num2str( sum(sum((imgsRef - imgsMov).^2)) / sizeImgs ) ]) ;
        end

% disp(['Deplacement ux : ', num2str(u(1)),' et deplacement uy : ', num2str(u(2)) , ' .' ]) ;
% disp(['Transformation Fxx = ', num2str(Am1(1,1)),', Fxy = ', num2str(Am1(1,2)) ]) ;
% disp(['               Fyx = ', num2str(Am1(2,1)) , ', Fyy = ', num2str(Am1(2,2)) , ' .' ]) 
%disp( [ 'Coefficient de Correlation (moyenne des carres des residus) : ' ,...
%    num2str( sum(sum((imgsRef - imgsMov).^2)) / sizeImgs ) ]) ;
%     % Timing
% disp(['ICGN: ',num2str(toc(startTime)*1000,'%.1f'),' ms']) ;

end