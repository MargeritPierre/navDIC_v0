function hd = captureInputs(hd)

    % Is there Inputs defined
        if isempty(hd.DAQInputs) ; return ; end
        if isempty(hd.DAQInputs.Inputs) ; return ; end
        
    % Acquire a single scan
        newData = inputSingleScan(hd.DAQInputs.Session).*[hd.DAQInputs.Inputs.Sensitivity] ;
        hd.InputData(hd.nFrames,1:length(newData)) = newData ;

end