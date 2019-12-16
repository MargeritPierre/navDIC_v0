function pts = my_ginput(PtsRef,ax1,ax2,CorrSize,siz)

pts = zeros(siz,2) ;

for i = 1 : siz
    ax1.XLim = [PtsRef(i,1)-2*CorrSize(1) , PtsRef(i,1)+2*CorrSize(1) ] ;
    ax1.YLim = [PtsRef(i,2)-2*CorrSize(2) , PtsRef(i,2)+2*CorrSize(2)] ;
    ax1.DataAspectRatio = [1,1,1] ;
    if i == 1 
        ax2.XLim = [PtsRef(i,1)-2*CorrSize(1) , PtsRef(i,1)+2*CorrSize(1)] ;
        ax2.YLim = [PtsRef(i,2)-2*CorrSize(2) , PtsRef(i,2)+2*CorrSize(2)] ;
    else
        ax2.XLim = [( PtsRef(i,1)-PtsRef(i-1,1)+pts(i-1,1) )-2*min(CorrSize) , ( PtsRef(i,1)-PtsRef(i-1,1)+pts(i-1,1) )+2*min(CorrSize)] ;
        ax2.YLim = [( PtsRef(i,2)-PtsRef(i-1,2)+pts(i-1,2) )-2*CorrSize(2) , ( PtsRef(i,2)-PtsRef(i-1,2)+pts(i-1,2) )+2*CorrSize(2)] ;
    end
    ax2.DataAspectRatio = [1,1,1] ;
    disp('Set the plot and press enter when ok') ;
    input('ready?') ;
    [ pts(i,1), pts(i,2) ] = ginput(1) ; 
end