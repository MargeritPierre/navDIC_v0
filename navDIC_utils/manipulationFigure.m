%% Rassemblement des resultats experimentaux 

figexp2 = gcf ;

PlotComp(3) = figexp2.Children(2).Children(1) ;
PlotComp(4) = figexp2.Children(2).Children(3) ;
PlotFil(3) = figexp2.Children(2).Children(2) ;
PlotFil(4) = figexp2.Children(2).Children(4) ;
%%
fig = figure ; 
plot( PlotFil(3).YData ,( PlotFil(3).XData + PlotFil(3).XData ) / 2 * 100)

%%
fig = figure ; 
plot( ( PlotFil(1).XData + PlotFil(2).XData ) / 2 * 100 - 0.022, PlotFil(1).YData )
plot( ( PlotComp(1).XData + PlotComp(2).XData ) / 2 * 100 + 0.01, PlotFil(1).YData )
legend({'\epsilon_{f} exp1','\epsilon_{c} exp 1'},'Interpreter','tex','location','Southeast') ;
xlabel('Deformation (%)')
ylabel('Force (N)')
title('Tension Stiffening ductal');

%%
fig = figure ; 
plot( ( PlotFil(1).XData + PlotFil(2).XData ) / 2 * 100 + 0.022, PlotFil(1).YData )
plot( ( PlotComp(1).XData + PlotComp(2).XData ) / 2 * 100 - 0.01, PlotFil(1).YData )
plot( ( PlotFil(3).XData + PlotFil(4).XData ) / 2 * 100 - 0.02, PlotFil(3).YData )
plot( ( PlotComp(3).XData + PlotComp(4).XData ) / 2 * 100 + 0.01, PlotFil(3).YData )
legend({'\epsilon_{f} exp1','\epsilon_{c} exp 1','\epsilon_{f} exp2','\epsilon_{c} exp 2'},'Interpreter','tex','location','Southeast') ;
xlabel('Deformation (%)')
ylabel('Force (N)')
title('Tension Stiffening ductal');
%%

for i = 1 : length(figexp.Children(2).Children)
    figexp.Children(2).Children(4).XData = figexp.Children(2).Children(4).XData*100 ;
end

%% Rassembler plusieurs plot
rep = 1 ;
i=0 ;
while rep
    disp('Selectionner un fichier figure matlab a assembler : ') ;
    [filename,path] = uigetfile('.fig','Selectionner un fichier figure matlab a assembler : ') ;
    fig = openfig([path,'/',filename]) ;
    for j = 1:length(fig.Children(2).Children)
        i = i+1 ;
        data(i).x = fig.Children(2).Children(j).XData ;
        data(i).y = fig.Children(2).Children(j).YData ;
        data(i).name = fig.Children(2).Children(j).DisplayName ;
    end
    rep = inputdlg('Voulez-vous fusionner une autre figure : ' ) ;
    rep = strcmpi(rep{1},'oui') ;
end
xlab = fig.Children(2).XLabel.String ;
ylab = fig.Children(2).YLabel.String ;

legtag = {} ;
% figure fusionn?e
figure ; 
xlabel(xlab,'Interpreter','tex');
ylabel(ylab,'Interpreter','tex');
for i = 1:length(data)
    plot(data(i).x,data(i).y) ;
    legtag{i} = data(i).name ;
end
legend(legtag,'Interpreter','tex','location','Southeast')

