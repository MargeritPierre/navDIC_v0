function  hd = useCameraCalibration(hd)
if length(hd.Cameras)<2 ; disp('you need to have 2 cameras'); return; end
if ~isfield( hd.Cameras, 'Properties' ); disp('Use hd = setCamerasProperties(hd,1)'); return; end


[filename,path] = uigetfile() ;

load( [path,filename] ) ;

% 
listSeed = {hd.Cameras.Name} ;
left_ID = listdlg('PromptString','Select the left Camera :',...
   'SelectionMode','Single',...
   'initialValue',1,...
   'ListString',listSeed) ;

if left_ID==1
    right_ID=2 ;
else
    right_ID=1 ;
end

% Change properties of cameras
hd.Cameras(left_ID).Properties.fpix = fc_left ;
hd.Cameras(right_ID).Properties.fpix = fc_right ;
hd.Cameras(left_ID).Properties.do = abs( T(1) ) ;
hd.Cameras(right_ID).Properties.do = abs( T(3) );
hd.Cameras(left_ID).Properties.pixObjratio = hd.Cameras(left_ID).Properties.fpix(1) / hd.Cameras(left_ID).Properties.do ;
hd.Cameras(right_ID).Properties.pixObjratio = hd.Cameras(right_ID).Properties.fpix(1) / hd.Cameras(right_ID).Properties.do ;


