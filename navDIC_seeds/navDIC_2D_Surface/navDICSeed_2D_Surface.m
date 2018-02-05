classdef navDICSeed_2D_Surface < navDICSeed
   
    properties
    end
    
    methods
        function obj = navDICSeed_2D_Surface(hd)
            % Initialize
                obj = obj@navDICSeed(hd,'single') ;
                obj.Class = 'navDICSeed_2D_Surface' ;
        end
    end
    
end
% 
% % Choose a camera
%         listCams = {hd.Cameras.Name} ;
%         [chosenID,valid] = listdlg('PromptString','Select Cameras to preview:',...
%                                     'SelectionMode','single',...
%                                     'initialValue',1,...
%                                     'ListString',listCams) ;
%         if ~valid ; return ; end
%         cam = hd.Cameras(chosenID) ;
%     
% 
%     % Take a picture
%         frame = step(cam.VidObj) ;
%         release(cam.VidObj) ;
% %% 
%     frame = zeros(480,640) ;
% 
%     % Set a ROI
%         [ROI,Shapes] = SetROI(frame,0) ;
%         Shapes ;
%         %%
%         clc
%     % Recursive Distance function
%         dF = {} ; % Individual shape dist function
%         dFrecurs = {} ; % Recursive dist function
%         pFix = [] ; % Fixed points
%         for s = 1:length(Shapes)
%             pos = Shapes(s).Position ;
%             switch Shapes(s).Class
%                 case 'imrect'
%                     dF{s} = @(p)drectangle(p,pos(1),pos(1)+pos(3),pos(2),pos(2)+pos(4)) ;
%                     pFix = [pFix ; ...
%                             pos(1) , pos(2) ;...
%                             pos(1) pos(2)+pos(4) ;...
%                             pos(1)+pos(3) pos(2) ;...
%                             pos(1)+pos(3) pos(2)+pos(4) ;...
%                             ] ;
%                 case 'imellipse'
%                     a = pos(3)/2 ;
%                     b = pos(4)/2 ;
%                     cx = pos(1)+a ;
%                     cy = pos(2)+b ;
%                     dF{s} = @(p)((p(:,1)-cx)./a).^2+((p(:,2)-cy)./b).^2-1 ;
%                 case 'impoly'
%                     dF{s} = @(p)dpoly(p,pos([1:end,1],:)) ;
%                     pFix = [pFix ; pos] ;
%             end
%             switch Shapes(s).Bool
%                 case '+'
%                     if s==1 
%                         dFrecurs{s} = @(p)dF{s}(p) ;
%                     else
%                         dFrecurs{s} = @(p)dunion(dFrecurs{s-1}(p),dF{s}(p)) ;
%                     end
%                 case '-'
%                     if s==1 
%                         dFrecurs{s} = @(p)dF{s}(p) ;
%                     else
%                         dFrecurs{s} = @(p)ddiff(dFrecurs{s-1}(p),dF{s}(p)) ;
%                     end
%             end
%         end
%     % Final dist. Function
%         fd = dFrecurs{end} ;
%     % Remove fixed points inside the domain
%         if ~isempty(pFix)
%             pFix = pFix(abs(fd(pFix))<=eps,:) ;
%         end
%     % BoundingBox
%         [j,i] = find(ROI) ;
%         margin = 1 ;
%         bboxROI = [min(i)-margin min(j)-margin ; max(i)+margin max(j)+margin] ;
%     % Compute the mesh
%         close all
%         figure('windowstyle','docked') ;
%         [Points,Tris] = distmesh2d(fd,@huniform,20,bboxROI,pFix);
%         