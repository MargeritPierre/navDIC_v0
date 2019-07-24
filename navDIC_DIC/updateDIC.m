function hd = updateDIC(hd)

    % Is there any DIC Seeds ?
        if isempty(hd.Seeds) ; return ; end
        
    % Is There any frame ?
        if hd.CurrentFrame==0 ; return ; end
        
    % Compute displacements
        for s = 1:length(hd.Seeds)
            if ~hd.Seeds(s).compDisp ; continue ; end
            hd.Seeds(s) = hd.Seeds(s).computeDisplacements(hd) ;
        end
        
    % Compute Strains
        for s = 1:length(hd.Seeds)
            if ~hd.Seeds(s).compDisp ; continue ; end
            if ~hd.Seeds(s).compStrains ; continue ; end
            hd.Seeds(s) = hd.Seeds(s).computeStrains(hd) ;
        end

end