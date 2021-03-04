%% WRITE A SERIES OF IMAGES as separated TXT files
global hd
IMG = {XX,YY} ; % hd.Images{3}
filenameTemplate = 'X%i.txt' ; 'Data/img_%04i.txt' ;
format = '%10.3f\t' ;

[nI,nJ,nC] = size(IMG{1}) ;
if nC>1 ; error('only grayscale images are supported') ; end
nIMG = numel(IMG) ;

wtbr = waitbar(0,'Exporting data to TXT file..') ;
for iii = 1:nIMG
    img = IMG{iii} ;
    filename = sprintf(filenameTemplate,iii) ;
    fid = fopen(filename,'w') ;
    fullFormat = [repmat(format,[1 nJ]) '\n'] ;
    fprintf(fid,fullFormat,img') ;
    fclose(fid) ;
    wtbr = waitbar(iii/nIMG,wtbr) ;
end
delete(wtbr) ;