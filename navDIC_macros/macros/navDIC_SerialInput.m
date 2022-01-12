classdef navDIC_SerialInput < navDIC_AbstractMacro
%NAVDIC_SerialInput Set an input as Serial Data

properties
    SerialPort internal.Serialport = internal.Serialport.empty
    SerialPortName char = 'none'
    BaudRate double = 115200
    LastMessage = "" ;
end

methods
    function this = navDIC_SerialInput()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
        this.SerialPort.configureCallback('off') ;
        delete(this.SerialPort);
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
    % Guess camID
    % Get parameters
        defInputs = { ...
                        'Input Name' , this.Name ...
                        ; strjoin(['Serial Port [none,' strjoin(serialportlist("all"),',') ']'],'') , this.SerialPortName ...
                        ; 'Baud Rate' , num2str(this.BaudRate) ...
                    } ;
        out = inputdlg(...
                    defInputs(:,1) ...
                    ,'Serial Input Properties' ...
                    ,1 ...
                    ,defInputs(:,2) ...
                    ) ;
    % Process output
        this.Name = out{1} ; 
        this.SerialPortName = out{2} ;
        this.BaudRate = str2double(out{3}) ;
    % Try openning the serial port
        if isvalid(this.SerialPort)
            this.SerialPort.configureCallback('off') ;
            delete(this.SerialPort) ;
        end
        if ~strcmp(this.SerialPortName,'none')
            try
                this.SerialPort = serialport(this.SerialPortName,this.BaudRate) ;
            catch 
                errordlg('Serial port failed to connect !','ERROR') ;
                return ;
            end
        end
    % Initialize the new input
        this.SerialPort.flush ;
        this.SerialPort.configureCallback('terminator',@this.serialCallback) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        data = str2num(this.LastMessage) ;
        data = data(:)' ;
        if isempty(hd.InputData) 
            hd.InputData = table() ; 
        end
        if ~ismember(this.Name,hd.InputData.Properties.VariableNames)
            hd.InputData.(this.Name) = NaN(hd.CurrentFrame,1) ;
        end
        if size(data,2)>size(hd.InputData.(this.Name),2)
            hd.InputData.(this.Name)(:,end+1:size(data,2)) = NaN ;
        end
        hd.InputData.(this.Name)(hd.CurrentFrame,1:size(data,2)) = data ;
    end
        
    function serialCallback(this,src,evt)
    % Callback function called on serial message received
        this.LastMessage = src.readline ;
        this.SerialPort.flush ;
    end
end

end

