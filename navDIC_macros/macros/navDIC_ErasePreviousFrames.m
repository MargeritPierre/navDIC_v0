classdef navDIC_ErasePreviousFrames < navDIC_AbstractMacro
%navDIC_ErasePreviousFrames Remove from hd.Images the old frames in
%order to reduce the RAM use.
% Only the most recent frame is kept in memory. 
properties
end

methods
    function this = navDIC_ErasePreviousFrames()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
    end
    
    function hd = setup(this,hd)
        %this.Name = inputdlg('Set the name of this CustomMacro:','Macro Name',1,{this.Name}) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        disp('Older frames removed from hd.Images') ;
        
    % Is there cameras ?
        if isempty(hd.Cameras) ; return ; end
        nCams = length(hd.Cameras) ;
        
    % Erase the latest image acquired and replace it by an empty cell
    for c = 1:nCams
        nFrames = length(hd.Images{c});
        if nFrames > 1
            if isa(hd.Images{c}{nFrames}, 'uint16')
                hd.Images{c}{nFrames - 1} = {uint16.empty(0)};
            elseif isa(hd.Images{c}{nFrames}, 'uint8')
                hd.Images{c}{nFrames - 1} = {uint8.empty(0)};
            else
                hd.Images{c}{nFrames - 1} =  {};
            end
        end
    end
    end
end

end

% Test : 
% >>> isempty(hd.Images{1,1}{1})
% ans = 1 

