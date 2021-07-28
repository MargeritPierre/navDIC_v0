global P8_Measures;

%% SENSITIVITY VALUES
forceS = 44.8;
dispS = 7.5;
LVDT_Supports_S = 2.5;
LVDT_MidSpan_S = 5.0;

%% DATA ACQ. FREQUENCIES
filterSpan = 1000; %DAQ
filterSpan_2 = 30; %OptiTracks

%% INITIALIZATION
P8_Measures.AcquiredData = [];
P8_Measures.AcquiredData.Data = [];
P8_Measures.AcquiredData.FilteredData = [];
P8_Measures.AcquiredData.InterpData = [];
P8_Measures.AcquiredData.InterpFilteredData = [];
P8_Measures.AcquiredData.Time = [];
P8_Measures.AcquiredData.StartTime = [];
P8_Measures.AcquiredOTPositions = [];
P8_Measures.AcquiredOTPositions.Data = [];
P8_Measures.AcquiredOTPositions.FilteredData = [];
P8_Measures.AcquiredOTPositions.InterpData = [];
P8_Measures.AcquiredOTPositions.InterpFilteredData = [];
P8_Measures.AcquiredOTPositions.InterpID = [];
P8_Measures.AcquiredOTPositions.Time = [];
P8_Measures.AcquiredOTPositions.ID = [];
P8_Measures.TimeLine = [];
P8_Measures.CommonTimeLine = [];
temp_Measures = [];
temp_Measures.Data = [];
temp_Measures.Time = [];

%% CUT THE DATAS WHEN ....time(end) IS DIFFERENT FOR SENSORS AND OT

[tMin_1, i1] = min([hd_1.AcquiredData.Time(end) - hd_1.AcquiredData.Time(1), ...
             hd_1.AcquiredOTPositions.Time(end) - hd_1.AcquiredOTPositions.Time(1)]);
         
[tMin_2, i2] = min([hd_2.AcquiredData.Time(end) - hd_2.AcquiredData.Time(1), ...
             hd_2.AcquiredOTPositions.Time(end) - hd_2.AcquiredOTPositions.Time(1)]);
 
[tMin_3, i3] = min([hd_3.AcquiredData.Time(end) - hd_3.AcquiredData.Time(1), ...
              hd_3.AcquiredOTPositions.Time(end) - hd_3.AcquiredOTPositions.Time(1)]);
        
if i1 == 1
    longestTimeStamps = hd_1.AcquiredOTPositions.Time ...
                        - hd_1.AcquiredOTPositions.Time(1);
    maxTimeIndex_1 = find((longestTimeStamps - tMin_1)>0);
    t1_1 = hd_1.AcquiredOTPositions.Time(maxTimeIndex_1(1) - 1) ... 
           - hd_1.AcquiredOTPositions.Time(1);
    t1_2 = hd_1.AcquiredOTPositions.Time(maxTimeIndex_1(1)) ... 
           - hd_1.AcquiredOTPositions.Time(1);
    hd_1.AcquiredOTPositions.Time = ...
        hd_1.AcquiredOTPositions.Time(1:maxTimeIndex_1(1));
    hd_1.AcquiredOTPositions.Time(end) = tMin_1 ...
                    + hd_1.AcquiredOTPositions.Time(1) ;
else
    longestTimeStamps = hd_1.AcquiredData.Time ...
                        - hd_1.AcquiredData.Time(1);
    maxTimeIndex_1 = find((longestTimeStamps - tMin_1)>0);
    t1_1 = hd_1.AcquiredData.Time(maxTimeIndex_1(1) - 1) ... 
           - hd_1.AcquiredData.Time(1);
    t1_2 = hd_1.AcquiredData.Time(maxTimeIndex_1(1)) ... 
           - hd_1.AcquiredData.Time(1);
    hd_1.AcquiredData.Time = ...
        hd_1.AcquiredData.Time(1:maxTimeIndex_1(1));
    hd_1.AcquiredData.Time(end) = tMin_1 ...
                    + hd_1.AcquiredData.Time(1);
end

if i2 == 1
    longestTimeStamps = hd_2.AcquiredOTPositions.Time ...
                        - hd_2.AcquiredOTPositions.Time(1);
    maxTimeIndex_2 = find((longestTimeStamps - tMin_2)>0);
    t2_1 = hd_2.AcquiredOTPositions.Time(maxTimeIndex_2(1) - 1) ... 
           - hd_2.AcquiredOTPositions.Time(1);
    t2_2 = hd_2.AcquiredOTPositions.Time(maxTimeIndex_2(1)) ... 
           - hd_2.AcquiredOTPositions.Time(1);
    hd_2.AcquiredOTPositions.Time = ...
        hd_2.AcquiredOTPositions.Time(1:maxTimeIndex_2(1));
    hd_2.AcquiredOTPositions.Time(end) = tMin_2 ...
                    + hd_2.AcquiredOTPositions.Time(1);
else
    longestTimeStamps = hd_2.AcquiredData.Time ... 
                        - hd_2.AcquiredData.Time(1);
    maxTimeIndex_2 = find((longestTimeStamps - tMin_2)>0);
    t2_1 = hd_2.AcquiredData.Time(maxTimeIndex_2(1) - 1) ... 
           - hd_2.AcquiredData.Time(1);
    t2_2 = hd_1.AcquiredData.Time(maxTimeIndex_2(1)) ... 
           - hd_2.AcquiredData.Time(1);
    hd_2.AcquiredData.Time = ...
        hd_2.AcquiredData.Time(1:maxTimeIndex_2(1));
    hd_2.AcquiredData.Time(end) = tMin_2 ...
                    + hd_2.AcquiredData.Time(1);
end

if i3 == 1
    longestTimeStamps = hd_3.AcquiredOTPositions.Time ...
                        - hd_3.AcquiredOTPositions.Time(1);
    maxTimeIndex_3 = find((longestTimeStamps - tMin_3)>0);
    t3_1 = hd_3.AcquiredOTPositions.Time(maxTimeIndex_3(1) - 1) ... 
           - hd_3.AcquiredOTPositions.Time(1);
    t3_2 = hd_3.AcquiredOTPositions.Time(maxTimeIndex_3(1)) ... 
           - hd_3.AcquiredOTPositions.Time(1);
    hd_3.AcquiredOTPositions.Time = ...
        hd_3.AcquiredOTPositions.Time(1:maxTimeIndex_3(1));
    hd_3.AcquiredOTPositions.Time(end) = tMin_3 ...
                    + hd_3.AcquiredOTPositions.Time(1);
else
    longestTimeStamps = hd_3.AcquiredData.Time ...
                        - hd_3.AcquiredData.Time(1);
    maxTimeIndex_3 = find((longestTimeStamps - tMin_3)>0);
    t3_1 = hd_3.AcquiredData.Time(maxTimeIndex_3(1) - 1) ... 
           - hd_3.AcquiredData.Time(1);
    t3_2 = hd_3.AcquiredData.Time(maxTimeIndex_3(1)) ... 
           - hd_3.AcquiredData.Time(1);
    hd_3.AcquiredData.Time = ...
        hd_3.AcquiredData.Time(1:maxTimeIndex_3(1));
    hd_3.AcquiredData.Time(end) = tMin_3 ...
                    + hd_3.AcquiredData.Time(1);
end

%% INTERPOLATE DATAS ACCORDING TO THE TIME CUTS

if i1 == 1
    hd_1.AcquiredOTPositions.Data = ...
                hd_1.AcquiredOTPositions.Data(:,:,1:size(hd_1.AcquiredOTPositions.Time,1));
            
    hd_1.AcquiredOTPositions.Data(:,:,end) = interpolateData(t1_1, t2_2, tMin_1, ...
                hd_1.AcquiredOTPositions.Data(:,:,end-1), ...
                hd_1.AcquiredOTPositions.Data(:,:,end));
            
    hd_1.AcquiredOTPositions.ID = ...
                hd_1.AcquiredOTPositions.ID(:,1:size(hd_1.AcquiredOTPositions.Time,1));
else
    hd_1.AcquiredData.Data = ...
                hd_1.AcquiredData.Data(1:size(hd_1.AcquiredData.Time,1),:);
            
    hd_1.AcquiredData.Data(end,:) = interpolateData(t1_1, t2_2, tMin_1, ...
                hd_1.AcquiredData.Data(end-1,:), ...
                hd_1.AcquiredData.Data(end,:));
end

if i2 == 1
    hd_2.AcquiredOTPositions.Data = ...
                hd_2.AcquiredOTPositions.Data(:,:,1:size(hd_2.AcquiredOTPositions.Time,1));
            
    hd_2.AcquiredOTPositions.Data(:,:,end) = interpolateData(t1_1, t2_2, tMin_1, ...
                hd_2.AcquiredOTPositions.Data(:,:,end-1), ...
                hd_2.AcquiredOTPositions.Data(:,:,end));
            
    hd_2.AcquiredOTPositions.ID = ...
                hd_2.AcquiredOTPositions.ID(:,1:size(hd_2.AcquiredOTPositions.Time,1));
else
    hd_2.AcquiredData.Data = ...
                hd_2.AcquiredData.Data(1:size(hd_2.AcquiredData.Time,1),:);
            
    hd_2.AcquiredData.Data(end,:) = interpolateData(t1_1, t2_2, tMin_1, ...
                hd_2.AcquiredData.Data(end-1,:), ...
                hd_2.AcquiredData.Data(end,:));
end

if i3 == 1
    hd_3.AcquiredOTPositions.Data = ...
                hd_3.AcquiredOTPositions.Data(:,:,1:size(hd_3.AcquiredOTPositions.Time,1));
            
    hd_3.AcquiredOTPositions.Data(:,:,end) = interpolateData(t1_1, t2_2, tMin_1, ...
                hd_3.AcquiredOTPositions.Data(:,:,end-1), ...
                hd_3.AcquiredOTPositions.Data(:,:,end));
            
    hd_3.AcquiredOTPositions.ID = ...
                hd_3.AcquiredOTPositions.ID(:,1:size(hd_3.AcquiredOTPositions.Time,1));
else
    hd_3.AcquiredData.Data = ...
                hd_3.AcquiredData.Data(1:size(hd_3.AcquiredData.Time,1),:);
            
    hd_3.AcquiredData.Data(end,:) = interpolateData(t1_1, t2_2, tMin_1, ...
                hd_3.AcquiredData.Data(end-1,:), ...
                hd_3.AcquiredData.Data(end,:));
end

%% ASSEMBLE DATAS FROM SEVERAL BATCHES

temp_Measures.Data = cat(1, hd_1.AcquiredData.Data, ...
                                       hd_2.AcquiredData.Data(2:end,:), ...
                                       hd_3.AcquiredData.Data(2:end,:));
                                   
temp_Measures.Time = cat(1, hd_1.AcquiredData.Time ...
                            - hd_1.AcquiredData.Time(1), ...
                            hd_2.AcquiredData.Time(2:end) ...
                            - hd_2.AcquiredData.Time(1) ...
                            + hd_1.AcquiredData.Time(end) ...
                            - hd_1.AcquiredData.Time(1), ...
                            hd_3.AcquiredData.Time(2:end) ...
                            - hd_3.AcquiredData.Time(1) ...
                            + hd_2.AcquiredData.Time(end) ...
                            - hd_2.AcquiredData.Time(1) ...
                            + hd_1.AcquiredData.Time(end) ...
                            - hd_1.AcquiredData.Time(1));
                     
P8_Measures.AcquiredData.StartTime = datetime(hd_1.AcquiredData.StartTime, ...
                                    'ConvertFrom','datenum', ...
                                    'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
                   
P8_Measures.AcquiredData.Time = P8_Measures.AcquiredData.StartTime ...
                                + seconds(temp_Measures.Time);
                            
P8_Measures.AcquiredOTPositions.Data = 1000.0 ...
                                * cat(3, hd_1.AcquiredOTPositions.Data, ...
                                hd_2.AcquiredOTPositions.Data(:,:,2:end), ...
                                hd_3.AcquiredOTPositions.Data(:,:,2:end));

P8_Measures.AcquiredOTPositions.ID = cat(2, hd_1.AcquiredOTPositions.ID, ...
                                           hd_2.AcquiredOTPositions.ID(:,2:end), ...
                                           hd_3.AcquiredOTPositions.ID(:,2:end));

P8_Measures.AcquiredOTPositions.Time = cat(1, hd_1.AcquiredOTPositions.Time ...
                                       - hd_1.AcquiredOTPositions.Time(1), ...
                                       hd_2.AcquiredOTPositions.Time(2:end) ...
                                       - hd_2.AcquiredOTPositions.Time(1) ...
                                       + hd_1.AcquiredOTPositions.Time(end) ...
                                       - hd_1.AcquiredOTPositions.Time(1), ...
                                       hd_3.AcquiredOTPositions.Time(2:end) ...
                                       - hd_3.AcquiredOTPositions.Time(1) ...
                                       + hd_2.AcquiredOTPositions.Time(end) ...
                                       - hd_2.AcquiredOTPositions.Time(1) ...
                                       + hd_1.AcquiredOTPositions.Time(end) ...
                                       - hd_1.AcquiredOTPositions.Time(1));
                                       
P8_Measures.AcquiredOTPositions.Time = P8_Measures.AcquiredData.StartTime ...
                                + seconds(P8_Measures.AcquiredOTPositions.Time);
                            
P8_Measures.TimeLine = cat(1, hd_1.TimeLine, hd_2.TimeLine, hd_3.TimeLine); 

P8_Measures.TimeLine = datetime(P8_Measures.TimeLine, ...
                                'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
                            
%% APPLY PROPER SENSITIVITY VALUES

forceMIN = min(temp_Measures.Data(:,1));
dispMAX = max(temp_Measures.Data(:,2));
disp_LVDT_1_MIN = min(temp_Measures.Data(:,3));
disp_LVDT_2_MIN = min(temp_Measures.Data(:,4));
disp_LVDT_3_MIN = min(temp_Measures.Data(:,5));

% forceSignal = (temp_Measures.Data(:,1) - forceMIN) * forceS;
% dispSignal = -(temp_Measures.Data(:,2) - dispMAX) * dispS;
% LVDT_1 = (temp_Measures.Data(:,3) - disp_LVDT_1_MIN) * LVDT_Supports_S;
% LVDT_2 = (temp_Measures.Data(:,4) - disp_LVDT_2_MIN) * LVDT_MidSpan_S;
% LVDT_3 = (temp_Measures.Data(:,5) - disp_LVDT_3_MIN) * LVDT_Supports_S;

forceSignal = temp_Measures.Data(:,1) * forceS;
dispSignal = -temp_Measures.Data(:,2) * dispS;
LVDT_1 = temp_Measures.Data(:,3) * LVDT_Supports_S;
LVDT_2 = temp_Measures.Data(:,4) * LVDT_MidSpan_S;
LVDT_3 = temp_Measures.Data(:,5) * LVDT_Supports_S;

P8_Measures.AcquiredData.Data = [forceSignal, dispSignal, LVDT_1, ...
                         LVDT_2, LVDT_3, temp_Measures.Data(:,6:end)];
                     
%% FILTER NON-INTERPOLATED SIGNALS

% TEST WITH GAUSSIAN FILTER
f1 = 3000;
w1 = gausswin(f1); %environ la fréquence d'acquisition
w1 = w1./sum(w1);

f2 = 30;
w2 = gausswin(f2); %environ la fréquence d'acquisition
w2 = w2./sum(w2);
% TEST WITH MOVING AVERAGE FILTERS
%movAvgFilter = 1/filterSpan * ones(filterSpan,1);
%delay_1 = floor(mean(grpdelay(movAvgFilter))); %Flooring is required for delay calibration
delay_1 = round(mean(grpdelay(w1)));

%movAvgFilter_2 = 1/filterSpan_2 * ones(filterSpan_2,1);
%delay_2 = floor(mean(grpdelay(movAvgFilter_2))); %Flooring is required for delay calibration
delay_2 = round(mean(grpdelay(w2)));

nbrOfSignals = size(P8_Measures.AcquiredData.Data, 2);
nbrOfMarkers = size(P8_Measures.AcquiredOTPositions.Data, 1);
            
P8_Measures.AcquiredData.FilteredData = filter(w1, 1, ...
               P8_Measures.AcquiredData.Data);

P8_Measures.AcquiredData.FilteredData(1:f1,:) = [];

P8_Measures.AcquiredData.Data = ...
                   P8_Measures.AcquiredData.Data(delay_1+1:end-delay_1,:);

P8_Measures.AcquiredData.Time = ...
                   P8_Measures.AcquiredData.Time(delay_1+1:end-delay_1,:);
            
P8_Measures.AcquiredOTPositions.FilteredData = filter(w2, 1, ...
                P8_Measures.AcquiredOTPositions.Data, [], 3);

P8_Measures.AcquiredOTPositions.FilteredData(:,:,1:f2) = [];

P8_Measures.AcquiredOTPositions.Data = ...
                   P8_Measures.AcquiredOTPositions.Data(:,:,delay_2+1:end-delay_2);

P8_Measures.AcquiredOTPositions.Time = ...
                   P8_Measures.AcquiredOTPositions.Time(delay_2+1:end-delay_2);
               
P8_Measures.AcquiredOTPositions.ID = ...
                   P8_Measures.AcquiredOTPositions.ID(:,delay_2+1:end-delay_2);

%% INTERPOLATE NON-FILTERED DATAS

%Using the interp1 function does not work on the timeline obtained directly
%because the timeline is not continuous. Plot the basic timelines to better
%understand.

Fr = 2; %Data frequency
P8_Measures.CommonTimeLine = P8_Measures.AcquiredData.Time(1):seconds(1/Fr):P8_Measures.AcquiredData.Time(end);

P8_Measures.AcquiredData.InterpData = ...
        interp1(P8_Measures.AcquiredData.Time(1:800:end).', ... %Supposing a frequency acquisition of 800Hz.
        P8_Measures.AcquiredData.Data(1:800:end,:), ...
        P8_Measures.CommonTimeLine);

R = reshape(P8_Measures.AcquiredOTPositions.Data(:,:,1:10:end), ...
        size(P8_Measures.AcquiredOTPositions.Data(:,:,1:10:end),1) ...
        * size(P8_Measures.AcquiredOTPositions.Data(:,:,1:10:end),2), []);

P8_Measures.AcquiredOTPositions.InterpData = ...
        interp1(P8_Measures.AcquiredOTPositions.Time(1:10:end), ... %Supposing a frequency acquisition of 30Hz
        R.', ...
        P8_Measures.CommonTimeLine);
    
P8_Measures.AcquiredOTPositions.InterpData = ...
        reshape(P8_Measures.AcquiredOTPositions.InterpData.', ...
        size(P8_Measures.AcquiredOTPositions.Data,1), ...
        size(P8_Measures.AcquiredOTPositions.Data,2),[]);

% %% FILTER INTERPOLATED SIGNALS
% 
% % TEST WITH GAUSSIAN FILTER
% winterp = gausswin(Fr); %environ la fréquence d'acquisition
% winterp = winterp./sum(winterp);
% 
% delay_interp = floor(mean(grpdelay(winterp)));
%             
% P8_Measures.AcquiredData.FilteredInterpData = filter(w1, 1, ...
%                 P8_Measures.AcquiredData.InterpData);
% 
% P8_Measures.AcquiredData.FilteredInterpData = ...
%                    P8_Measures.AcquiredData.FilteredInterpData((Fr+1):end,:);
%                
% P8_Measures.AcquiredData.FI_Time = P8_Measures.AcquiredData.Time((Fr+1):end); %Shrink the time vector to the same length as the data ones. 
% 
% P8_Measures.AcquiredOTPositions.FilteredInterpData = filter(Fr, 1, ...
%                 cat(3, P8_Measures.AcquiredOTPositions.InterpData, ...
%                 zeros(nbrOfMarkers,3,Fr)), [], 3);
%             
% P8_Measures.AcquiredOTPositions.FilteredInterpData = ...
%                    P8_Measures.AcquiredOTPositions.FilteredInterpData(:,:,(Fr+1):end);

%% "NORMALIZE" UNFILTERED NON-INTERPOLATED SIGNALS WITH FIRST DATA

% forceMIN = P8_Measures.AcquiredData.Data(1,1);
% dispMAX = P8_Measures.AcquiredData.Data(1,2);
% disp_LVDT_1_MIN = P8_Measures.AcquiredData.Data(1,3);
% disp_LVDT_2_MIN = P8_Measures.AcquiredData.Data(1,4);
% disp_LVDT_3_MIN = P8_Measures.AcquiredData.Data(1,5);
% 
% forceSignal = P8_Measures.AcquiredData.Data(:,1) - forceMIN;
% dispSignal = P8_Measures.AcquiredData.Data(:,2) - dispMAX;
% LVDT_1 = P8_Measures.AcquiredData.Data(:,3) - disp_LVDT_1_MIN;
% LVDT_2 = P8_Measures.AcquiredData.Data(:,4) - disp_LVDT_2_MIN;
% LVDT_3 = P8_Measures.AcquiredData.Data(:,5) - disp_LVDT_3_MIN; 
% 
% P8_Measures.AcquiredData.Data = [forceSignal, dispSignal, LVDT_1, ...
%                          LVDT_2, LVDT_3, P8_Measures.AcquiredData.Data(:,6:end)];
%             
% P8_Measures.AcquiredOTPositions.Data = P8_Measures.AcquiredOTPositions.Data ...
%                       - P8_Measures.AcquiredOTPositions.Data(:,:,1);

%% "NORMALIZE" FILTERED SIGNALS WITH FIRST DATA

% forceMIN = P8_Measures.AcquiredData.FilteredData(1,1);
% dispMAX = P8_Measures.AcquiredData.FilteredData(1,2);
% disp_LVDT_1_MIN = P8_Measures.AcquiredData.FilteredData(1,3);
% disp_LVDT_2_MIN = P8_Measures.AcquiredData.FilteredData(1,4);
% disp_LVDT_3_MIN = P8_Measures.AcquiredData.FilteredData(1,5);
% 
% forceSignal = P8_Measures.AcquiredData.FilteredData(:,1) - forceMIN;
% dispSignal = P8_Measures.AcquiredData.FilteredData(:,2) - dispMAX;
% LVDT_1 = P8_Measures.AcquiredData.FilteredData(:,3) - disp_LVDT_1_MIN;
% LVDT_2 = P8_Measures.AcquiredData.FilteredData(:,4) - disp_LVDT_2_MIN;
% LVDT_3 = P8_Measures.AcquiredData.FilteredData(:,5) - disp_LVDT_3_MIN; 
% 
% P8_Measures.AcquiredData.FilteredData = [forceSignal, dispSignal, LVDT_1, ...
%                          LVDT_2, LVDT_3, P8_Measures.AcquiredData.FilteredData(:,6:end)];
%             
% P8_Measures.AcquiredOTPositions.FilteredData = P8_Measures.AcquiredOTPositions.FilteredData ...
%                       - P8_Measures.AcquiredOTPositions.FilteredData(:,:,1);
                  
%% "NORMALIZE" FILTERED INTERPOLATED SIGNALS WITH FIRST DATA

forceMIN = P8_Measures.AcquiredData.InterpData(1,1);
dispMAX = P8_Measures.AcquiredData.InterpData(1,2);
disp_LVDT_1_MIN = P8_Measures.AcquiredData.InterpData(1,3);
disp_LVDT_2_MIN = P8_Measures.AcquiredData.InterpData(1,4);
disp_LVDT_3_MIN = P8_Measures.AcquiredData.InterpData(1,5);

forceSignal = P8_Measures.AcquiredData.InterpData(:,1) - forceMIN;
dispSignal = P8_Measures.AcquiredData.InterpData(:,2) - dispMAX;
LVDT_1 = P8_Measures.AcquiredData.InterpData(:,3) - disp_LVDT_1_MIN;
LVDT_2 = P8_Measures.AcquiredData.InterpData(:,4) - disp_LVDT_2_MIN;
LVDT_3 = P8_Measures.AcquiredData.InterpData(:,5) - disp_LVDT_3_MIN; 

P8_Measures.AcquiredData.InterpData = [forceSignal, dispSignal, LVDT_1, ...
                         LVDT_2, LVDT_3, P8_Measures.AcquiredData.InterpData(:,6:end)];
            
P8_Measures.AcquiredOTPositions.InterpData = P8_Measures.AcquiredOTPositions.InterpData ...
                      - P8_Measures.AcquiredOTPositions.InterpData(:,:,1);
                  
%% INTERPOLATE FILTERED DATAS

%Using the interp1 function does not work on the timeline obtained directly
%because the timeline is not continuous. Plot the basic timelines to better
%understand.

% Fr = 2; %Data frequency
% P8_Measures.CommonTimeLine = P8_Measures.AcquiredData.Time(1):seconds(1/Fr):P8_Measures.AcquiredData.Time(end); 

P8_Measures.AcquiredData.InterpFilteredData = ...
        interp1(P8_Measures.AcquiredData.Time(1:800:end).', ... %Supposing a frequency acquisition of 800Hz
        P8_Measures.AcquiredData.FilteredData(1:800:end,:), ...
        P8_Measures.CommonTimeLine);

R = reshape(P8_Measures.AcquiredOTPositions.FilteredData(:,:,1:10:end), ...
        size(P8_Measures.AcquiredOTPositions.FilteredData(:,:,1:10:end),1) ...
        * size(P8_Measures.AcquiredOTPositions.FilteredData(:,:,1:10:end),2), []);

P8_Measures.AcquiredOTPositions.InterpFilteredData = ...
        interp1(P8_Measures.AcquiredOTPositions.Time(1:10:end), ... %Supposing a frequency acquisition of 30Hz
        R.', ...
        P8_Measures.CommonTimeLine);
    
P8_Measures.AcquiredOTPositions.InterpFilteredData = ...
        reshape(P8_Measures.AcquiredOTPositions.InterpFilteredData.', ...
        size(P8_Measures.AcquiredOTPositions.FilteredData,1), ...
        size(P8_Measures.AcquiredOTPositions.FilteredData,2),[]);

P8_Measures.AcquiredOTPositions.ID = double(P8_Measures.AcquiredOTPositions.ID);

Rid = P8_Measures.AcquiredOTPositions.ID(:,1:10:end);
    
P8_Measures.AcquiredOTPositions.InterpID = ...
        interp1(P8_Measures.AcquiredOTPositions.Time(1:10:end), ... %Supposing a frequency acquisition of 30Hz
        Rid.', ...
        P8_Measures.CommonTimeLine);
    
P8_Measures.AcquiredOTPositions.InterpID = P8_Measures.AcquiredOTPositions.InterpID.';
    
%% USEFUL FUNCTION

% Interpolating function
function [interpData] = interpolateData(t1, t2, tD, D1, D2)
    interpData = D1 + (tD - t1)/(t2 - t1) * D2;
end
