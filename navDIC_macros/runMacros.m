function hd = runMacros(hd,callback)
% Run macros on a navDIC event
% "callback" allows to choose the kind of function to execute
if nargin<2 ; callback = 'run' ; end

    % Is there Macros defined ?
        if isempty(hd.Macros) ; return ; end
        if ~isa(hd.Macros,'navDIC_AbstractMacro') ; return ; end
        
    % Run the macro's corresponding callback
        for m = 1:numel(hd.Macros)
            hd = hd.Macros(m).(callback)(hd) ;
        end

end