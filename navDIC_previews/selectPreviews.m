function [IDs,valid] = selectPreviews(hd,SelectionMode)
    % Init
        IDs = [] ;
        valid = false ;
    % Are there any previews ?
        if isempty(hd.Previews) ; return ; end
    % Is there only one preview ?
        if length(hd.Previews)==1 
            IDs  = 1 ; 
            valid = true ; 
            return ;
        end
    % Otherwise, choose a preview
        listPreviews = [hd.Previews{:}] ;
        listPreviews = [listPreviews.fig] ;
        listPreviews = {listPreviews.Name} ;
        [IDs,valid] = listdlg('PromptString','Select a Preview :',...
                                    'SelectionMode',SelectionMode,...
                                    'initialValue',1,...
                                    'ListString',listPreviews) ;
end