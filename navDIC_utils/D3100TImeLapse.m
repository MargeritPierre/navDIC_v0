%% Identification
% Science Dude 1990
% December 10, 2020
%
%% Code
% Take a photo with a Nikon D3100 and transfer to local PC, for making
% timelapse movies!!!

%% Clean up
clc
close all
clear
drawnow

%% Parameters
% Pause after all the other tasks (i.e., the main timelapse parameter)
pause_tasks = 0;

% Number of pictures to take
N_pictures = 5;

% COM port
com_port = 5;

% BAUD rate
baud_rate = 9600;

% Data bits
data_bits = 8;

% Camera directory
camera_dir = 'Computer\D3100\Removable storage\DCIM\100TEST_';

% Current directory
current_dir = pwd;

% The filename
temp = clock;

% Filename for saving
f_name = ['Timelapse_Record_', num2str(temp(1), '%04d'), '_', num2str(temp(2), '%02d'), '_', num2str(temp(3), '%02d'), '_', ...
          num2str(temp(4), '%02d'), '_', num2str(temp(5), '%02d'), '_', num2str(round(temp(6)), '%02d')];

%% Open the COM port
% The serial port
s = serial(['COM' num2str(com_port)], 'BaudRate', baud_rate, 'DataBits', data_bits, 'Terminator', 'CR');
fopen(s);
pause(1);
% Save the response from the microcontroller (opening com port resets the
% controller)
scan_count = 1;
tmp = {};
while s.BytesAvailable > 1
    tmp{scan_count} = fscanf(s);
    scan_count = scan_count + 1;
end

%% Active X to be able to copy the pictures from the camera to the local PC
% Create active x server
h = actxserver('WScript.Shell');

%% Main loop
% Cell array to hold the COM port communications
com_record = cell(N_pictures, 3);

% The pause to let windows react
pause_ui = 0.15;

% Pause for the snapshot
pause_snapshot = 1.5;

for ii = 1 : N_pictures
    %% Send the COM commands
    pause(pause_ui);
    % Read the register for reference
    fprintf(s, 'r37');
    pause(pause_ui);
    % Set PORTB to 1, then 3, to take a picture
    fprintf(s, 'w37=1');
    pause(pause_ui);
    fprintf(s, 'w37=3');
    % Note the time the picture was requested
    com_record{ii, 1} = clock;
    pause(pause_ui);
    fprintf(s, 'w37=0');
    
    % Let the picture get taken
    pause(pause_snapshot);
    
    % Save the response from the microcontroller
    scan_count = 1;
    tmp = {};
    while s.BytesAvailable > 1
        tmp{scan_count} = fscanf(s);
        scan_count = scan_count + 1;
    end
    % Save for later, for review if needed
    com_record{ii, 2} = tmp;
    
    %% Go get the photo from camera
        
    % Counters
    picture_ready_count = 0;
    picture_ready = 0;
    
    % In case the picture fails to take, place some text in the clipboard
    clipboard('copy', 'NOTAFILENAME');
        
    while picture_ready == 0
        % Raise the camera explorer view
        temp = h.AppActivate(camera_dir);
    
        if temp ~= 1
            error('Could not raise camera directory');
        end
    
        % Work around to get focus to the picture in the camera
        h.SendKeys('%{TAB}');
        pause(pause_ui);
        h.SendKeys('%{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
    
        % Send the "down" key
        h.SendKeys('{DOWN}');                   
        pause(pause_ui);
        % Copy the filename of the new picture to the clipboard
        h.SendKeys('{F2}');
        pause(pause_ui);
        h.SendKeys('^c');
        % Let the filename get picked up by the clipboard
        pause(pause_ui);
        % Get the filename of the picture from the clipboard
        pic_filename = clipboard('paste');
    
        if isempty(pic_filename) || strcmp(pic_filename, 'NOTAFILENAME')
            picture_ready = 0;
            picture_ready_count = picture_ready_count + 1;
            pause(0.25);
        else            
            % Save the filename in the record
            com_record{ii, 3} = pic_filename;
            
            % Try to copy the file
            pause(pause_ui);
            h.SendKeys('{esc}');
            pause(pause_ui);
            h.SendKeys('{esc}');
            pause(pause_ui);
            
            % "Cut" the file from the camera
            h.SendKeys('^x');
            pause(pause_ui);

            h.SendKeys('^x');
            pause(pause_ui);
            
            % Check MATLAB's view of the clipboard - if a file is "cut"
            % then the clipboard response to paste should be empty
            temp = clipboard('paste');
            
            if isempty(temp)
                picture_ready = 1;
            end
        end
        
        if picture_ready_count > 30
            error('Picture taking error');
        end
    end    
    
    % Make sure the file got copied
    wait_count = 0;
    wait_for_file = 1;
    while wait_for_file == 1
        % Raise the local folder
        temp = h.AppActivate(current_dir);
    
        if temp ~= 1
            error('Could not raise local directory');
        end

        % Work around to get focus to the directory
        h.SendKeys('%{TAB}');
        pause(pause_ui);
        h.SendKeys('%{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{TAB}');
        pause(pause_ui);
        h.SendKeys('{DOWN}');
        pause(pause_ui);

        % "Paste" the picture from the camera
        h.SendKeys('^v');
        % Let the file get moved
        pause(0.25);
   
        temp = dir([pic_filename '.*']);
        
        if isempty(temp)
            wait_for_file = 1;
            wait_count = wait_count + 1;
            pause(0.25);
            
        else
            wait_for_file = 0;
        end
        
        if wait_count > 20
            error('File did not copy');
        end
    end    
        
    %% Final pause
    pause(pause_tasks);

    % Save the record of the timelapse      
    save(f_name, 'com_record');  
    
    disp(['Photo ' num2str(ii) ' of ' num2str(N_pictures)]);
    
end
    
%% Clean up
% Close the serial port
fclose(s);