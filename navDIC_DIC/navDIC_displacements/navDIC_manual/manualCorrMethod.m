function MovingPoints = manualCorrMethod(PtsMov,PtsRef,imgMov,imgRef,CorrSize)
    
    % Convertir image en 8 bits en float 
        if isa(imgMov(1,1),'uint8')
            imgMov = double(imgMov)/255 ;
        end
        if isa(imgRef(1,1),'uint8')
            imgRef = double(imgRef)/255 ;
        end

    % PARAMETERS
        fig = figure ; 
        ax1 = subplot(1,2,1);
        ax1.YDir = 'reverse' ; 
        image(repmat(imgRef,[1,1,3])) ; 
        plot(PtsRef(:,1), PtsRef(:,2) ,'sq','LineWidth',1,'MarkerSize',10)
        siz = size( PtsRef , 1) ; 
        txt = {} ;
        for i = 1 : siz ; txt{i} = ['\leftarrow ', num2str(i)] ; end
        text( PtsRef(:,1), PtsRef(:,2),txt,'interpreter','tex' )
        
        ax2 = subplot(1,2,2);
        ax2.YDir = 'reverse' ; 
        image(repmat(imgMov,[1,1,3])) ; 
        MovingPoints = my_ginput(PtsRef,ax1,ax2,CorrSize,siz) ;
        close(fig)
end