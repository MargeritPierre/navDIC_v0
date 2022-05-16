

    % CLEAN WORKSPACE
        clc
        global hd
        clearvars -except hd
        
    % PARAMETERS
        seedNumber = 1 ;
        frames = [21:27] ;
        Neig = 3 ;
        normalize = true ;
        
    % GET THE DATA
        Nodes = hd.Seeds(seedNumber).Points ;
        Elems = hd.Seeds(seedNumber).Triangles ;
        U = hd.Seeds(seedNumber).Displacements(:,:,frames) ;
        
    % GRADIENT OPERATOR (U->Strains)
        Xt = reshape(Nodes(Elems',1),3,[])' ; Yt = reshape(Nodes(Elems',2),3,[])' ;
        areas = 1/2*(Xt(:,2).*Yt(:,3) - Xt(:,3).*Yt(:,2) - Xt(:,1).*Yt(:,3) + Xt(:,3).*Yt(:,1) + Xt(:,1).*Yt(:,2) - Xt(:,2).*Yt(:,1)) ;
        aa = cross(Xt,Yt,2) ;
        bb = circshift(Yt,2,2) - circshift(Yt,1,2) ;
        cc = circshift(Xt,1,2) - circshift(Xt,2,2) ;
        Dx = sparse(repmat((1:size(Elems,1))',[1 3]),Elems,1./2./areas(:).*bb,size(Elems,1),size(Nodes,1)) ;
        Dy = sparse(repmat((1:size(Elems,1))',[1 3]),Elems,1./2./areas(:).*cc,size(Elems,1),size(Nodes,1)) ;
        E = zeros(size(Elems,1),3,length(frames)) ;
        E(:,1,:) = Dx*reshape(U(:,1,:),size(Nodes,1),[]) ;
        E(:,2,:) = Dy*reshape(U(:,2,:),size(Nodes,1),[]) ;
        E(:,3,:) = Dy*reshape(U(:,1,:),size(Nodes,1),[]) + Dx*reshape(U(:,2,:),size(Nodes,1),[]) ;
        
    % PERFORM THE EVD
        DATA = U ;
        S = reshape(DATA,[],length(frames)) ;
        if normalize 
            NORM = sqrt(sum(S.^2,2)) ;
        else
            NORM = ones(size(S,1),1) ;
        end
        S = S./NORM ;
        if size(S,1)<size(S,2)
            Css = S*S' ;
            [V,d] = eigs(Css,Neig,'lm') ;
            d = sqrt(diag(d)) ;
            T = S'*V*diag(1./d) ;
        else
            Css = S'*S ;
            [T,d] = eigs(Css,Neig,'lm') ;
            d = sqrt(diag(d)) ;
            V = S*T*diag(1./d) ;
        end
        V = V.*NORM ;
        V = reshape(V,size(DATA,1),size(DATA,2),Neig) ;
        MDL = abs(diff(d))./d(1:end-1) ;
        
        if 0
            E = zeros(size(Elems,1),3,Neig) ;
            E(:,1,:) = Dx*reshape(V(:,1,:),size(Nodes,1),[]) ;
            E(:,2,:) = Dy*reshape(V(:,2,:),size(Nodes,1),[]) ;
            E(:,3,:) = Dy*reshape(V(:,1,:),size(Nodes,1),[]) + Dx*reshape(V(:,2,:),size(Nodes,1),[]) ;
            V = E ;
        end
        
    % PLOT THE RESULT
    
        colorV = @(i) sqrt(sum(V(:,:,i).^2,2)) ;
        mode = 1 ;
        
        clf ;
        if size(V,1) == size(Nodes,1)
            FaceColor = 'interp' ;
        else
            FaceColor = 'flat' ;
        end
        axD = axes('nextplot','add',) ;
            %plot(cumsum(d)/sum(d),':ok') ; axD.XLim = [1 Neig] ; axD.YLim = [0.9 1] ; axD.YScale = 'log' ; grid on
            %plot(d,':ok') ; axD.XLim = [1 Neig] ; axD.YScale = 'log' ; grid on
            plot(T(:,mode)*d(mode),'k') ;
        ax = gobjects(0) ;
        meshV = gobjects(0) ;
        for c = 1:size(V,2) 
            ax(c) = axes('nextplot','add','position',[1 1 0 0] + 0.5*[-1 0 1 0] + 0.3*[0 -c*1.05 0 1]) ;
                %patch('Vertices',Nodes,'Faces',Elems,'facecolor','none','edgecolor','k','edgealpha',0.5) ;
                meshV(c) = patch('Vertices',Nodes,'Faces',Elems,'facecolor',FaceColor,'facealpha',0.5,'edgecolor','none','edgealpha',0.5) ;
                meshV(c).FaceVertexCData = V(:,c,mode) ;
                if size(V,1) == size(Nodes,1)
                    meshV(c).Vertices(:,3) = V(:,c,mode) ;
                else
                    %meshV(c).Vertices(:,3) = V(:,c,mode) ;
                end
        end
        %axis(ax,'equal')
        axis(ax,'tight')
        axis(ax,'off')
        
        
        
        
        