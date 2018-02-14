classdef navDIC_Source < handle & matlab.mixin.Heterogeneous
    
    properties
        % Class, Name and validity
            Name = '' ;
            Class = 'BasicSource' ;
            isValid = true ;
        % Kepp all Acquired data ? 
            KeepData = true ;
        % Acquired Data
            firstFrame = true ; % Only one frame has been pushed
            Data = [] ;
    end
    
    
    methods
        
        % CONSTRUCTOR
            function obj = navDIC_Source()
            end
        
        % GET DATA AT SPECIFIC FRAMES
            function data = getData(obj,frames)
                if obj.KeepData && ~obj.firstFrame
                    % Sub referencing
                        s.type = '()' ;
                        s.subs = cellstr(repmat(':',[ndims(obj.Data)-1 ,1]))' ;
                        s.subs{end+1} = frames ;
                        data = subsref(obj.Data,s) ;
                else
                    if length(frames)>1
                        data = [] ;
                    else
                        data = obj.Data ;
                    end
                end
            end
            
        % GET LAST ACQUIRED DATA
            function data = getLastData(obj)
                if obj.firstFrame
                    data = obj.Data ;
                else
                    data = getData(obj,size(obj.Data,ndims(obj.Data))) ;
                end
            end
            
        % GET ALL ACQUIRED DATA
            function data = getAllData(obj)
                data = obj.Data ;
            end
        
        % PUSH DATA (ADD DATA AT THE END)
            function pushData(obj,data)
                if obj.KeepData
                    obj.Data = cat(ndims(data)+1,obj.Data,data) ;
                else
                    obj.Data = data ;
                end
                if size(obj.Data,ndims(data)+1)>1
                    obj.firstFrame = false ;
                end
            end
        
        % CLEAR DATA
            function clearData(obj)
                obj.Data = [] ;
                obj.firstFrame = true ;
            end
            
    end
    
end