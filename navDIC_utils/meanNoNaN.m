function out = meanNoNaN(data,dim)
    isNotNAN = double(~isnan(data)) ;
    data(~isNotNAN) = 0 ;
    out = sum(data,dim)./sum(isNotNAN,dim) ;