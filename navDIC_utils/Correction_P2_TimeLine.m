%% hd_2
ADCorrected = [];
ADCorrected.Time = [];
ADCorrected.Data = [];
j = 0;
n = 0;
m = 0;
for i=1:length(hd_2.AcquiredData.Time)
    if(hd_2.AcquiredData.Time(i) == 10*j)
        if(n==2)
            n=0;
            j = j+1;
        else
            n = n+1;
        end
    end
    if (hd_2.AcquiredData.Time(i) >= 10*j && n == 2)
        ADCorrected.Time = [ADCorrected.Time hd_2.AcquiredData.Time(i)];
        ADCorrected.Data = cat(1,ADCorrected.Data,hd_2.AcquiredData.Data(i,:));
    end
end  

hd_2.AcquiredData.Time = ADCorrected.Time;
hd_2.AcquiredData.Data = ADCorrected.Data;
plot(ADCorrected.Time)
save('D:\00_THESE\03_PRODUCTION\03_OPTIMISATION\02_POUTRE EUROCODE\04_Analyse Structurale\Essais\Flexion\22032021_Test_P2_PoutreOptimisée\hd_2_CorrectedTime.mat','hd_2');

%% hd_3
ADCorrected = [];
ADCorrected.Time = [];
ADCorrected.Data = [];
j = 0;
n = 0;
m = 0;
for i=1:length(hd_3.AcquiredData.Time)
    if(hd_3.AcquiredData.Time(i) == 10*j)
        if(n==2)
            n=0;
            j = j+1;
        else
            n = n+1;
        end
    end
    if (hd_3.AcquiredData.Time(i) >= 10*j && n == 2)
        ADCorrected.Time = [ADCorrected.Time hd_3.AcquiredData.Time(i)];
        ADCorrected.Data = cat(1,ADCorrected.Data,hd_3.AcquiredData.Data(i,:));
    end
end  

hd_3.AcquiredData.Time = ADCorrected.Time;
hd_3.AcquiredData.Data = ADCorrected.Data;
plot(ADCorrected.Time)
save('D:\00_THESE\03_PRODUCTION\03_OPTIMISATION\02_POUTRE EUROCODE\04_Analyse Structurale\Essais\Flexion\22032021_Test_P2_PoutreOptimisée\hd_3_CorrectedTime.mat','hd_3');

%% hd_4
ADCorrected = [];
ADCorrected.Time = [];
ADCorrected.Data = [];
j = 0;
n = 0;
m = 0;
for i=1:length(hd_4.AcquiredData.Time)
    if(hd_4.AcquiredData.Time(i) == 10*j)
        if(n==2)
            n=0;
            j = j+1;
        else
            n = n+1;
        end
    end
    if (hd_4.AcquiredData.Time(i) >= 10*j && n == 2)
        ADCorrected.Time = [ADCorrected.Time hd_4.AcquiredData.Time(i)];
        ADCorrected.Data = cat(1,ADCorrected.Data,hd_4.AcquiredData.Data(i,:));
    end
end  

hd_4.AcquiredData.Time = ADCorrected.Time;
hd_4.AcquiredData.Data = ADCorrected.Data;
plot(ADCorrected.Time)
save('D:\00_THESE\03_PRODUCTION\03_OPTIMISATION\02_POUTRE EUROCODE\04_Analyse Structurale\Essais\Flexion\22032021_Test_P2_PoutreOptimisée\hd_4_CorrectedTime.mat','hd_4');

%% hd_5
ADCorrected = [];
ADCorrected.Time = [];
ADCorrected.Data = [];
j = 0;
n = 0;
m = 0;
for i=1:length(hd_5.AcquiredData.Time)
    if(hd_5.AcquiredData.Time(i) == 10*j)
        if(n==2)
            n=0;
            j = j+1;
        else
            n = n+1;
        end
    end
    if (hd_5.AcquiredData.Time(i) >= 10*j && n == 2)
        ADCorrected.Time = [ADCorrected.Time hd_5.AcquiredData.Time(i)];
        ADCorrected.Data = cat(1,ADCorrected.Data,hd_5.AcquiredData.Data(i,:));
    end
end  

hd_5.AcquiredData.Time = ADCorrected.Time;
hd_5.AcquiredData.Data = ADCorrected.Data;
plot(ADCorrected.Time)
save('D:\00_THESE\03_PRODUCTION\03_OPTIMISATION\02_POUTRE EUROCODE\04_Analyse Structurale\Essais\Flexion\22032021_Test_P2_PoutreOptimisée\hd_5_CorrectedTime.mat','hd_5');

%% hd_6
ADCorrected = [];
ADCorrected.Time = [];
ADCorrected.Data = [];
j = 0;
n = 0;
m = 0;
for i=1:length(hd_6.AcquiredData.Time)
    if(hd_6.AcquiredData.Time(i) == 10*j)
        if(n==2)
            n=0;
            j = j+1;
        else
            n = n+1;
        end
    end
    if (hd_6.AcquiredData.Time(i) >= 10*j && n == 2)
        ADCorrected.Time = [ADCorrected.Time hd_6.AcquiredData.Time(i)];
        ADCorrected.Data = cat(1,ADCorrected.Data,hd_6.AcquiredData.Data(i,:));
    end
end  

hd_6.AcquiredData.Time = ADCorrected.Time;
hd_6.AcquiredData.Data = ADCorrected.Data;
plot(ADCorrected.Time)
save('D:\00_THESE\03_PRODUCTION\03_OPTIMISATION\02_POUTRE EUROCODE\04_Analyse Structurale\Essais\Flexion\22032021_Test_P2_PoutreOptimisée\hd_6_CorrectedTime.mat','hd_6');

%% hd_7
ADCorrected = [];
ADCorrected.Time = [];
ADCorrected.Data = [];
j = 0;
n = 0;
m = 0;
for i=1:length(hd_7.AcquiredData.Time)
    if(hd_7.AcquiredData.Time(i) == 10*j)
        if(n==2)
            n=0;
            j = j+1;
        else
            n = n+1;
        end
    end
    if (hd_7.AcquiredData.Time(i) >= 10*j && n == 2)
        ADCorrected.Time = [ADCorrected.Time hd_7.AcquiredData.Time(i)];
        ADCorrected.Data = cat(1,ADCorrected.Data,hd_7.AcquiredData.Data(i,:));
    end
end  

hd_7.AcquiredData.Time = ADCorrected.Time;
hd_7.AcquiredData.Data = ADCorrected.Data;
plot(ADCorrected.Time)
save('D:\00_THESE\03_PRODUCTION\03_OPTIMISATION\02_POUTRE EUROCODE\04_Analyse Structurale\Essais\Flexion\22032021_Test_P2_PoutreOptimisée\hd_7_CorrectedTime.mat','hd_7');
