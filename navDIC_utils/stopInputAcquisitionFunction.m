function stopInputAcquisitionFunction()

% FUNCTION EXECUTED AT THE END OF THE TIMER.
    global hd   
    global plotVar  
    stop(hd.DAQInputs.Session);
    hd.OT.disable(0);
    stop(plotVar.plotTimer); 

end