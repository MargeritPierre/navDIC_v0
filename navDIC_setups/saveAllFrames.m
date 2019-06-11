function hd = saveAllFrames(hd) 

    for f = 1:hd.nFrames
        hd = saveCurrentSetup(hd,f) ;
    end