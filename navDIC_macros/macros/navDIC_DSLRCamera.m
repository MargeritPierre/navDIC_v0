classdef navDIC_DSLRCamera < navDIC_VirtualCamera
%NAVDIC_DSLRCamera Virtual camera with DSLR camera

properties
    SerialPort %internal.Serialport = internal.Serialport.empty
    SerialPortName char = 'none'
    BaudRate double = 115200
    Folder char = ''
    TimeOut double = 10
end

methods
    function this = navDIC_DSLRCamera()
    % Class constructor
    end
    
    function delete(this)
    % Class destructor
        if ~isempty(this.SerialPort) && isvalid(this.SerialPort)
            fclose(this.SerialPort) ;
            delete(this.SerialPort) ;
        end
    end
    
    function hd = setup(this,hd)
    % Setup the macro (change parameters, etc)
    % Guess camID
        if isempty(this.CameraID) ; this.CameraID = numel(hd.Cameras)+1 ; end
    % Get parameters
        defInputs = { ...
                        'Name' , this.Name ...
                        ; 'ID' , num2str(this.CameraID) ...
                        ; strjoin(['Serial Port [none,' strjoin(seriallist("all"),',') ']'],'') , this.SerialPortName ...
                        ; 'Baud Rate' , num2str(this.BaudRate) ...
                        %; 'Image Folder' , this.Folder ...
                        %; 'Time out' , num2str(this.TimeOut) ...
                    } ;
        out = inputdlg(...
                    defInputs(:,1) ...
                    ,'DSLR Camera Properties' ...
                    ,1 ...
                    ,defInputs(:,2) ...
                    ) ;
    % Process output
        this.Name = out{1} ; 
        this.SerialPortName = out{3} ;
        this.BaudRate = str2double(out{4}) ;
        %this.Folder = out{5} ;
        %this.TimeOut = str2double(out{6}) ;
    % Try openning the serial port
        if ~isempty(this.SerialPort) && isvalid(this.SerialPort)
            fclose(this.SerialPort) ;
            delete(this.SerialPort) ;
        end
        if ~strcmp(this.SerialPortName,'none')
            try
                this.SerialPort = serial(this.SerialPortName,'BaudRate',this.BaudRate) ;
                fopen(this.SerialPort) ;
            catch 
                errordlg('Serial port failed to connect !','ERROR') ;
                return ;
            end
        end
    % Try getting a new frame
%         try
%             img = this.getFrame() ;
%         catch 
%             errordlg('Failed to get a frame !','ERROR') ;
%             return ;
%         end
%         if isempty(img) ; errordlg('No frame acquired !','ERROR') ; return ; end
    % Set the navDIC virtual cam
        %this.ImageSize = size(img) ;
        %hd = this.setCam(hd,str2double(out{2})) ;
    end

    function hd = run(this,hd)
    % Function executed when a new frame is added to navDIC
        %hd = run@navDIC_VirtualCamera(this,hd) ;
        %if isempty(this.CameraID) ; return ; end
        img = this.getFrame ;
        %if isempty(img) ; return ; end
        %hd.Images{this.CameraID}{hd.CurrentFrame} = img ;
    end
    
    function img = getFrame(this)
    % Get a frame
        img = [] ;
    % Check that everything is OK
        if isempty(this.SerialPort) || ~isvalid(this.SerialPort) ; return ; end
        %if isempty(this.Folder) || ~isfolder(this.Folder) ; return ; end
    % Get the folder state
        %files0 = dir(this.Folder) ;
    % Send the shoot flag
        fprintf(this.SerialPort,"shoot") ;
    % Wait for the folder to change
        %t = tic ;
        %while toc(t)<this.TimeOut
        %end
    end
end

end

