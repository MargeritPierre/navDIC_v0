function hd = runMacros(hd)

    % Is there Macros defined ?
        if isempty(hd.Macros) ; return ; end
        
    % Run the macro
        for m = 1:numel(hd.Macros)
            hd = hd.Macros(m).run(hd) ;
        end

end