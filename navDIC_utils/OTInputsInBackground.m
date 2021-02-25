global hd
optTrack = hd.OT;

if exist ('lst','var') && isvalid(lst) ; delete (lst) ; end
optTrack.addlistener(1, 'Position');
optTrack.enable(1);
%pause(3)
%optTrack.disable(0);

% global hd
% hd.AcquiredOTPositions = [];
% hd.AcquiredOTPositions.Positions = [];
% hd.AcquiredOTPositions.Time = [];
% 
% function onOTDataAvailable(src, event)
%     global hd
%     event;
%     hd.onOTDataAvailable.Positions = [hd.onOTDataAvailable.Positions ; event];
% end