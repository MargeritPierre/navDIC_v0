%% PRE ALLOCATE FRAME FOR REAL-TIME NAVDIC
    global hd
    nFr = 300 ; 

    for c = 1:numel(hd.Cameras) 
        cam = hd.Cameras(c).VidObj ;
        roi = cam.ROIPosition ;
        im = getsnapshot(cam) ;
        hd.Images{c} = num2cell(zeros([roi([4 3]) 1 nFr],class(im)),1:3) ; 
    end