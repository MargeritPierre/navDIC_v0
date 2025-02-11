%% CREATE A navDIC SESSION WITH EVERYTHING READY
% Before, a reference frame needs to be present in the watched folder

% Initialization
cd 'E:\In-situ\MATLAB' % folder that contains the startup & navDIC .m files
startup

% Choose the working directory
defaultWD = 'E:\Partage_MEB\Temp_Movefile_image' ;
wd = uigetdir(defaultWD,'CHOOSE THE WORKING DIRECTORY') ;
if wd==0 ; return ; end

% Restart navDIC
navDIC exit ;
cd(wd) ;
navDIC ;
global hd ;

% Watch folder macro
disp("SETUP THE WATCH FOLDER MACRO")
hd.Macros = navDIC_WatchFolder() ;
hd.Macros.WatchFolder = wd ; % watch the working directory by default
hd.Macros.TimeOut = 60 ;
hd.Macros.Delay = 0.5 ;
hd = hd.Macros.setup(hd,1) ; % setup so that it is the first camera
hd.script.ui.update() ; % update the navDIC toolbar

% Create the first seed
disp("CREATE THE FIRST POINT SEED")
hd.CurrentFrame = 1 ; % set the reference frame to the first frame
hd.Seeds = navDICSeed_2D_DistMesh(hd) ;
hd.Seeds.Name = 'RefPoint' ;

% Setup the FFT (local) DIC
disp("SETUP THE LOCAL DIC")
hd.Macros(2) = navDIC_FFTDICMacro() ;
hd.Macros(2).Seed = hd.Seeds(1) ;
hd.Macros(2).UsePreviousVelocity = 0 ;
hd.Macros(2).setup(hd) ;
hd.Macros(2).CorrSize = round(0.5*autoCorrSize(hd.Macros(2))) ; % auto corrsize choice

% GlobalDIC mesh(es)
disp("CREATE THE MESH SEED(S)")
hd.Seeds(2:end) = [] ;
while true
    try ; hd.Seeds(end+1) = navDICSeed_2D_DistMesh(hd) ;
    catch ; break ;
    end
    hd.Seeds(end).Name = char("Mesh"+string(numel(hd.Seeds)-1)) ;
end

% GlobalDIC Macros
disp("CREATE THE GLOBAL DIC MACROS")
hd.Macros(3:end) = [] ;
for ss = 1:numel(hd.Seeds)-1
% Projection macro
    hd.Macros(end+1) = navDIC_ProjectMesh() ;
    hd.Macros(end).MasterSeed = hd.Seeds(ss) ;
    hd.Macros(end).SlaveSeed = hd.Seeds(ss+1) ;
    hd.Macros(end).initProjector() ;
% Global DIC macro
    hd.Macros(end+1) = navDIC_GlobalDICMacro() ;
    hd.Macros(end).Seed = hd.Seeds(ss+1) ;
    hd.Macros(end).UsePreviousVelocity = 0 ;
    hd.Macros(end).MaxDisp = 0.01 ; % convergence criterion !
    hd.Macros(end).setup(hd) ;
end

% Update the toolbar
hd.script.ui.update() ; % update the navDIC toolbar

% Re-run macros (will perform the DIC(s) on existing images)
hd.script.ui.clickOn('reRunMacros') ;

% Seed Preview
hd.script.ui.clickOn('seedPreview') ;

% Plot Preview
hd.script.ui.clickOn('plotPreview') ;