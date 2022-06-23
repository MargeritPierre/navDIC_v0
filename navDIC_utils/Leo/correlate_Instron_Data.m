% This scripts helps recorrelating data given by the instron table w/
% something from the DIC measurement, if any has failed. For example, u
% could recorrelate the start of displacement identified thru DIC to the
% plate displacement, and then get back the force values.

% First : compute the interplate displacement with navDIC
% 


%% DONT EXECUTE ALL THE SCRIPT, just run the sections accordingly

%% Save the datas.

results.pixelMmRatio(test) = rapport;
results.ysPeakMin(test) = ysPeakMin;
results.ysPeakMax(test) = ysPeakMax;
results.ysPlateauMin(test) = ysPlateauMin;
results.ysPlateauMax(test) = ysPlateauMax;
results.ysPeakDiamMax(test) = peakMaxDiameterPixel;
results.ysPeakDiamMin(test) = peakMinDiameterPixel;

results.ysPlateauDiamMin(test) = plateauMinDiameterPixel;
results.ysPlateauDiamMax(test) = plateauMaxDiameterPixel;

results.peakFrame(test) = peakFrame;
results.plateauFrame(test) = plateauFrame;

%% 1) Load the DIC data. Initialize everything for a new test treatment
global hd
FillMissingImages;
test = 1;

path_to_adamel = 'C:\Users\Leo\Documents\MATLAB\2202-02-SqueezeTest2\donnees_adamel\1.txt'
table_machine = readtable(path_to_adamel);

%% Set variables to compute static (peak) and dynamic (plateau) yield stress.
peakFrame = 29;
plateauFrame = 32;


%% Set the scaling factor (u might use the commented part to recompute it)
% Check image scale
image = hd.Images{1,1}{1,peakFrame};
DefineScale(image);

%% Check diameter for plateauframe
image = hd.Images{1,1}{1,plateauFrame};
DefineScale(image);

%%
plateauMaxDiameterPixel = 2838.85;
plateauMinDiameterPixel = 2750.33;
peakMaxDiameterPixel = 2739.02;
peakMinDiameterPixel = 2739.02;

pixelMMRatio = 148/3703.63; % About 0.0400
%% Load Datas.
% Load the instron data table.

% Go to navDIC and compute a displacement (u2).


%% Correlate data
% Compute the offset (Offs variable)


%% Offset, resample machine data and save as hd variable (inputdata)

% Load the DIC data u want to recorrelate (e.g. a distmesh between two
% plates)
dic_disp_pixels = squeeze(hd.Seeds.DataFields.u2);
dic_disp_pixels = dic_disp_pixels(4,:);
frequency_dic = 2;
%dic_disp = dic_disp(1,:);
% identify table framerate.
timestep = table_machine.Temps_s_(4)-table_machine.Temps_s_(3);
frequency_machine = 1/timestep;
% Knowing the framerate, resample the DIC values to display as time series
% if necessary.


% 
disp_dic_mm = -(dic_disp_pixels)*pixelMMRatio;


disp_machine_mm = table_machine.Allongement_mm_-table_machine.Allongement_mm_(1);
force_machine_n = table_machine.Force_N_-table_machine.Force_N_(1);

time_dic = (1:length(disp_dic_mm))/frequency_dic;
time_machine = (1:length(disp_machine_mm))/frequency_machine;

%% Compute offset
offs = MoveCrv(disp_dic_mm,time_dic,disp_machine_mm,time_machine)


%% Resample instron data as DIC data.
% Offs is to be found before with the MoveCrv function.
interpMachineDisp = griddedInterpolant(time_machine+offs,disp_machine_mm);
interpMachineForce = griddedInterpolant(time_machine+offs,force_machine_n);
interpMachineForce.Method = 'spline';
interpMachineForce.ExtrapolationMethod = 'nearest';
resampledMachineDisp = interpMachineDisp(time_dic);
resampledForce = interpMachineForce(time_dic);

hd.InputData(:,2) = resampledForce';
figure
plot(time_machine+offs,force_machine_n);
plot(time_dic,resampledForce);



%% Compute peak yield stress

ysPeakMax = YS(peakFrame,peakMaxDiameterPixel,pixelMMRatio);
ysPeakMin = YS(peakFrame,peakMinDiameterPixel,pixelMMRatio);

disp(strcat('peakYsMax : ', num2str(ysPeakMax),' Pa'));
disp(strcat('peakysmin : ', num2str(ysPeakMin),' Pa'));

ysPlateauMax = YS(plateauFrame,plateauMaxDiameterPixel,pixelMMRatio);
ysPlateauMin = YS(plateauFrame,plateauMinDiameterPixel,pixelMMRatio);

disp(strcat('plateauYsMax : ', num2str(ysPlateauMax),' Pa'));
disp(strcat('plateauysmin : ', num2str(ysPlateauMin),' Pa'));


%%

% % Calculation data for Sample1
% pixeldiameterAtPlateau1 = 2804.81;
% % pixelDiameterAtPeak1 = 2668.61;
% pixeldiamPeak2 = 2655.00;
% pixeldiamPlateau2 = 2891.09;
% diameterAtPlateau = pixeldiameterAtPlateau * pixelMMRatio;
% d = pixeldiamPlateau2 * pixelMMRatio;
% % diameterAtPeak1 = pixelDiameterAtPeak1*pixelMMRatio;
% % d = diameterAtPeak1;
% transitionFrameIndex = 82;
% %% COMPUTE YIELD STRESS
% transitionRadius = (d/2);
% areaTransition = pi*transitionRadius^2;
% stressAtTransition = (hd.InputData(transitionFrameIndex,2) / areaTransition)*1e+6; % Pa
% disp(strcat('stress : ', num2str(stressAtTransition),' Pa'));





%plot(time_machine+10,disp_machine_mm);

% figure

% plot(time_machine+offs,disp_machine);
% plot(time_dic,disp_dic1);
% plot(time_dic,disp_dic2);
% plot(time_dic,disp_dic3);
% plot(time_dic,disp_dic4);
% plot(time_dic,disp_dic5);
% legend('machine','4,:',"1,:","mean 4 and 1", "median","mean");

% figure;
% title('raw displacement (pixel) vs time')
% plot(time_dic,-dic_disp)

% figure;
% plot((1:length(disp_machine))/frequency_machine,disp_machine);
% plot((1:length(disp_dic))/frequency_dic,disp_dic);
% legend('machine','dic')
% Display both graphs as normalized values, and offset the table values to
% match the DIC values.

% Create the function to dynamically offset machine value.

% preview recorrelated force value.

% Now find the closest index for both sets.

% Cross-identify force value.

function offs = MoveCrv(dic_disp,dic_time,machine_disp,machine_time)

hold on
fig = figure;

set(fig, 'KeyPressFcn', @keypress);

hax = axes('Parent',fig);
offs = 0;

plot(dic_time,dic_disp, 'DisplayName', 'Machine','Parent',hax);

x = machine_time;
y = machine_disp;
h = plot(x,y, ':', 'DisplayName', 'DIC','Parent',hax);

  redraw();
    
waitfor(fig);
%             e = SaveData(e);
            disp(strcat('recorded refit. Coefficient is : ',num2str(offs,1)))
    

function keypress(~, evnt)
        switch lower(evnt.Key)
            case 'rightarrow'
                offs = offs+0.055;
            case 'leftarrow'
                offs = offs-0.055;
            otherwise
                return
        end
        % Always do a redraw
        redraw();
    end
        
    function redraw()
        set(h, 'XData', x+offs, ...
               'YData', y);
        drawnow
    end

%     function e = SaveData(e)
%     e.dic.strain.refit = x+offs;
%     e.dic.strain.refitCoef = offs;
%     end


end



%% Subfunctions

function ys = YS(frame,diameter_px,ratio)
global hd
peakDiameterMm = diameter_px * ratio;
dpe = peakDiameterMm;
transitionRadius = (dpe/2);
areaTransition = pi*transitionRadius^2;
ys = (hd.InputData(frame,2) / areaTransition)*1e+6; % Pa
end

function scale = DefineScale(image)
figure;
a = imshow(image);
% start drawing
h = imdistline;
api = iptgetapi(h);
dist = api.getDistance();
fprintf('The length of the segment is: %0.2f pixels \n', dist)
scale = dist;
end