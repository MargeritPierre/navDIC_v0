function [C,lvl] = tricontours(t,p,f,lvl)

% Deal with given levels
    if nargin<4
        lvl = 10 ;
    end
    if numel(lvl)==1
        lvl = linspace(1/(lvl+1),1-1/(lvl+1),lvl)*range(f(:)) + min(f(:)) ;
    end
    
% Edges
    e = cat(3,t,circshift(t,1,2)) ;
% Function values at edges nodes
    f1 = f(e(:,:,1)) ;
    f2 = f(e(:,:,2)) ;
% Edges nodes positions
    x1 = reshape(p(e(:,:,1),1),[],3) ;
    y1 = reshape(p(e(:,:,1),2),[],3) ;
    x2 = reshape(p(e(:,:,2),1),[],3) ;
    y2 = reshape(p(e(:,:,2),2),[],3) ;
    
% Contour lines
    C = [] ;
    for ll = 1:numel(lvl)
        L = lvl(ll) ;
        % Parameter of intersection (t € [0,1])
            t = (L-f1)./(f2-f1) ;
        % Points of intersection
            xc = x1.*(1-t) + x2.*t ;
            yc = y1.*(1-t) + y2.*t ;
        % Indices of edges crossing the z=L plane (valid parameters)
            iE = t<1 & t>0 & ~isnan(t) ;
        % Indices triangles crossing the z=L plane (two edges cross)
            iT = sum(iE,2)==2 ;
        % Make cross lines
            [iT,iE] = find(iE & iT) ;
            [iT,ind] = sort(iT) ; iE = iE(ind) ;
            iP = sub2ind(size(t),iT,iE) ;
            X = reshape(xc(iP),2,[]) ;
            Y = reshape(yc(iP),2,[]) ;
        % Contour lines
            nans = X(1,:)*NaN ;
            X = [X ; nans] ;
            Y = [Y ; nans] ;
            Z = X*0 + L ;
            C = [C cat(3,X,Y,Z)] ;
    end

