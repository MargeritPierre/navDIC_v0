function stopInputAcquisitionFunction()

% FUNCTION EXECUTED AT THE END OF ThdHE TIMER.
    global hd   
    global plotVar 
    
    hd.AcquiredData.EndTime = [];
    hd.AcquiredOTPositions.EndTime = [];
    
    hd.AcquiredData.EndTime = now;
    stop(hd.DAQInputs.DataAcquisition);
    flush(hd.DAQInputs.DataAcquisition);
    %stop(hd.DAQInputs.Session);
    hd.AcquiredOTPositions.EndTime = now;
    hd.OT.disable(0);
    %stop(plotVar.plotTimer);
end