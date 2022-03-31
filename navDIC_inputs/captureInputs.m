function hd = captureInputs(hd)

    % Is there Inputs defined
        if isempty(hd.DAQInputs) ; return ; end
        if isempty(hd.DAQInputs.Inputs) ; return ; end
        
    % Acquire a single scan
        try
            newData = inputSingleScan(hd.DAQInputs.Session) ;
        catch
            if ~isfield(hd.DAQInputs,'Background')
                hd.DAQInputs.Background = struct() ;
                hd.DAQInputs.Background.Listener = addlistener(hd.DAQInputs.Session,'DataAvailable',@dataAvailableCallback) ;
                hd.DAQInputs.Background.LastData = NaN(1,numel(hd.DAQInputs.Inputs)) ;
                hd.DAQInputs.Session.IsContinuous = true ;
                hd.DAQInputs.Session.startBackground() ;
            end
            newData = hd.DAQInputs.Background.LastData ;
        end
        newData = newData.*[hd.DAQInputs.Inputs.Sensitivity] ;
        hd.InputData(hd.nFrames,1:length(newData)) = newData ;

end

function dataAvailableCallback(src,evt)
    global hd
    hd.DAQInputs.Background.LastData = evt.Data(end,:) ;
end