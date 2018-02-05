%% IMAGE ACQUISITION AND CORRELATION UTILITARY

clear all
close all
clc
        
%% SET SESSION FOR EXTERNAL INPUTS (FORCE SIGNAL, ETC..)
    clc
    % Create session if needed
        daqreset ; clear('session') ;
        if ~exist('session','var') || ~isvalid(session)
            session = daq.createSession('ni') ;
        end
    % Get available devices
        dev = daq.getDevices() ;
    % Allow the choice of a device
        nD = length(dev);
        switch nD
            case 0 % No Devices
                display([char(10),'NO AVAILABLE DEVICE']) ; 
                return
            case 1 % Only One Device
                num = 1 ;
            otherwise % A choice has to be made
            display([char(10), 'AVAILABLE DEVICES :'])
            for d = 1:nD
                display(['   ',num2str(d),' : ',dev(d).Vendor.FullName,' ',dev(d).Model]) ;
            end
            num = input('Choose a device : ') ;
        end
        dev = dev(num) ;
    % Allow to choose the device subsystem
        nS = length(dev.Subsystems);
        switch nS
            case 0
                display([char(10),'NO SUBSYSTEMS']) ; 
                return
            case 1
                num = 1 ;
            otherwise % A choice has to be made
            display([char(10),dev.Vendor.FullName,' ',dev.Model, ', subsystems :']) ;
            for s = 1:nS
                display(['   ',num2str(s),' : ',dev.Subsystems(s).SubsystemType]) ;
            end
            num = input('Choose a subsystem : ') ;
        end
        sub = dev.Subsystems(num) ;
    % Allow to choose an input
        nI = length(sub.ChannelNames);
        switch nI
            case 0
                display([char(10),'NO INPUTS AVAILABLE']) ; 
                return
            case 1
                num = 1 ;
            otherwise % A choice has to be made
            display([char(10),sub.SubsystemType,' available inputs :']) ;
            for i = 1:nI
                display(['   ',num2str(i-1),' : ',sub.ChannelNames{i}]) ;
            end
            num = input('Choose an input : ') ;
        end
        in = num+1 ;
    % Add the input to the session and NAME IT
        in = addAnalogInputChannel(session,dev.ID,sub.ChannelNames{in},'Voltage') ;
        in.Name = input('Specify a name for the input : ','s') ;
    % USER-DEFINED PROPERTIES
        in.Range = [-4,4];
        Sensivity = 4953 ; % N/V
        U0 = inputSingleScan(session)*Sensivity ;
        Units = 'N' ;
    % Display the session
        display(session) ;
    
%% DISPLAY THE INPUT
    fig = figure('windowstyle','docked') ;
    fig.Position(3:4) = [1 1]*400 ;
    DATA = inputSingleScan(session)*Sensivity-U0 ;
    pl = plot(0,0,'k') ;
    ttl = title('bla','fontweight','bold') ;
    t = tic ;
    while ishandle(fig)
        data = inputSingleScan(session);
        DATA(end+1,:) = data*Sensivity-U0 ;
        try
            pl.YData = [pl.YData , DATA(end,:)] ;
            pl.XData = [pl.XData , toc(t)] ; (1:size(DATA,1)-1)'*ones(1,size(DATA,2)) ;
            ttl.String = [in.Name,' : ',num2str(DATA(end,:),4),' ',Units] ;
        end
        drawnow ;
        pause(0.05) ;
    end 



%% CONNECT TO A CAMERA AND CHOOSE THE IMAGE FORMAT
% info on dcam camera on canal 1

info = imaqhwinfo('dcam',1) ;
formats = info.SupportedFormats' ;

    nF = length(info.SupportedFormats);
    display([char(10), 'AVAILABLE FORMATSS :'])
    for f = 1:nF
        display(['   ',num2str(f),' : ',formats{f}])
    end
    num = input('Choose a format : ');

    while or(num>nF, num<0)
        num = input('entrer le numero de la resolution choisie: ');
        disp(num);
    end
    reso = info.SupportedFormats(num);
    vidobj = imaq.VideoDevice('dcam', 1, reso{1});
    
    vidobj.ReturnedColorSpace = 'grayscale';
    vidobj.ReturnedDataType = 'uint8';
    vidobj.DeviceProperties.ShutterMode = 'manual';
    vidobj.DeviceProperties.GainMode = 'manual';
    %vidobj.ROI = [1 336 1600 324];
    
    
%% PREVIEW CAMERA OUTPUT WITH ALLOWED ZOOM (FOR LENS AND BRIGTNESS TUNING)
    obj = preview(vidobj) ;
    obj.Parent.YDir = 'reverse' ;
    zoom(obj.Parent,'on') ;

    
%% SET DATA ACQUISITION FOLDER
    [file,path] = uiputfile('*','SELECTIONNER REPERTOIRE DE SAUVEGARDE ET NOM IMAGES','img.tif') ;
    if file ==0 ; return ; end
    [~,file,ext] = fileparts(file) ;
  
    
%% SET REGION OF INTEREST
    % Take a picture as refImg for ROI setting
        roiImg = step(vidobj);
    % Open the GUI for ROI setting
        ROI = SetROI(roiImg) ;
    % Save the ROI as image
        fileNameROI = [path,file,'_ROI',ext] ;
        imwrite(ROI,fileNameROI);
    % Get the boundingBox of ROI
        [i,j] = find(ROI) ;
        bboxROI = [min(i) min(j) max(i)-min(i) max(j)-min(j)] ;

%% ACQUIRE IMAGES AT A FIXED FRAMERATE
    clc
    close all
        
    % Parameters
        FrameRate = 1 ; % Hz
        saveImages = true ;
        % Cropping
            CropWithROI = true ;
            marginX = 0 ; % pixels margin to include in the cropping
            marginY = 0 ; % pixels margin to include in the cropping
            
    % Cropping zone
        frame = step(vidobj) ;
        if CropWithROI
            imin = max(1,bboxROI(1)-marginY) ;
            jmin = max(1,bboxROI(2)-marginX) ;
            imax = min(size(frame,1),bboxROI(1)+bboxROI(3)+marginY) ;
            jmax = min(size(frame,2),bboxROI(2)+bboxROI(4)+marginX) ;
        else 
            imin = 1 ; imax = size(frame,1) ;
            jmin = 1 ; jmax = size(frame,2) ;
        end
        iCrop = imin:imax ;
        jCrop = jmin:jmax ;
        frame = frame(iCrop,jCrop,:) ;
        
    % Initialize figure
        btn_Height = 40 ;
        btn_Width = 80 ;
        btn_Margin = 10 ;
        fig = figure ;
            fig.Position([4,3]) = size(frame)+(2*btn_Margin+btn_Height)*[1 0] ;
            ax = axes('position',[0 0 1 1]) ;
            ax.Units = 'pixels' ;
            ax.Position = ax.Position + [0 2*btn_Margin+btn_Height 0 -2*btn_Margin-btn_Height] ;
            ax.YDir = 'reverse' ;
            im = imagesc(frame) ;
            colormap(gray)
            axis equal
            axis tight
            axis off
        btn = uicontrol('style','togglebutton','position',[btn_Margin btn_Margin btn_Width btn_Height]) ;
            btn.String = 'START' ;
            btn.FontWeight = 'bold' ;
            btn.FontSize = 15 ;
        info = uicontrol('style','text','position',[btn.Position(1)+btn.Position(3)+btn_Margin btn_Margin fig.Position(3)*10 btn_Margin+btn_Height/2]) ;
            info.BackgroundColor = 'w' ;
            info.HorizontalAlignment = 'left' ;
            info.FontSize = 10 ;
            info.FontWeight = 'bold' ;
            info.String = '' ;
        drawnow ;
            
            
    % BEGGIN THE ACQUISITION
        % Wait until the START button is Clicked
            while ishandle(fig) && ~btn.Value
                drawnow ;
            end
            btn.String = 'STOP' ;
        % Start Acquisition
            im_id = 0 ;
            time = [] ;
            data = [] ;
            t = tic ;
            while ishandle(fig) && btn.Value
                % Increase image ID
                    im_id = im_id+1 ;
                % Get the time 
                    time(end+1) = toc(t) ;
                % Get a frame and crop it
                    frame = step(vidobj) ;
                    frame = frame(iCrop,jCrop,:) ;
                % Get input data if needed
                    if exist('session','var') && isvalid(session)
                        data(end+1,:) = inputSingleScan(session).*Sensivity-U0 ;
                    end
                % Save the image
                    if saveImages 
                        fileName = [path,file,'_',num2str(im_id),ext] ;
                        imwrite(frame, fileName) ;
                    end
                % Update figure 
                    try % if the figure is closed, no errors
                        im.CData = frame ;
                        info.String = [num2str(im_id),' images'] ;
                        if exist('session','var') && isvalid(session) % prompt input values
                            info.String = [info.String,', ',in.Name,' : ',num2str(data(end,:),4),' ',Units] ;
                        end
                    end
                    drawnow ;
                % Wait until next Img
                    while toc(t)*FrameRate<im_id
                        drawnow ;
                    end
            end
        % When Acquisition Stopped
            % Disable Button
                btn.Enable = 'off' ;
                btn.String = 'END' ;
                csvData = [] ;
            % Save TimeLine
                if saveImages
                    save([path,file,'_time.mat'],'time')
                    csvData.Headers{1} = 'time' ;
                    csvData.Data(:,1) = time(:) ;
                end
            % Save input data if needed
                if saveImages && exist('session','var') && isvalid(session)
                    save([path,file,'_inputs.mat'],'data') ;
                    csvData.Headers{end+1} = [in.Name,' (',Units,')'] ;
                    csvData.Data(:,end+1) = data ;
                end
            % Save the CROPPED ROI as image
                if saveImages
                    fileNameROI = [path,file,'_ROI',ext] ;
                    imwrite(ROI(iCrop,jCrop,:),fileNameROI);
                end
            % Export as CSV
                if saveImages
                    csvData.FileName = [path,file,'_CSV.csv'] ;
                    fid = fopen(csvData.FileName, 'w') ;
                    % Headers
                        if length(csvData.Headers)>1
                            fprintf(fid, '%s;', csvData.Headers{1:end-1}) ;
                        end
                        fprintf(fid, '%s;\n', csvData.Headers{end}) ;
%                     % Data
                            %dlmwrite(csvData.FileName, csvData.Data+eps, '-append') ;
                            %csvwrite(csvData.FileName,csvData.Data,0,0) ;
                        for i = 1:length(time)
                            if length(csvData.Data(i,:))>1
                                fprintf(fid, '%s;', csvData.Data(i,1:end-1)+eps) ;
                            end
                            fprintf(fid, ['%s;',char(10)], csvData.Data(i,end)+eps) ;
                        end
                    % Close file
                    fclose(fid) ;
                    %
                end
                
                
