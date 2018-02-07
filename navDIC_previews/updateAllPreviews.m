function hd = updateAllPreviews(hd)

    % Is there previews ?
        if isempty(hd.Previews) ; return ; end
        
    % Update previews
        validPrev = [] ;
        for p = 1:length(hd.Previews)
            hd.Previews{p} = hd.Previews{p}.updatePreview(hd) ;
            validPrev(p) = hd.Previews{p}.isValid ;
        end
        validPrev = logical(validPrev) ;
        
    % Clear all invalid previews
        %if all(validPrev) ; return ; end
        if ~any(validPrev) 
            hd.Previews = {} ;
            return ; 
        end
        hd.Previews = hd.Previews(validPrev) ;

end