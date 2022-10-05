% D=xlsread('synergyfinder.xlsx');
D=xlsread('SynergyFinderPlus.csv');

%% jcnkn
A_P_theo=D(1:64,6);
A_T_theo=D(65:128,6);
P_T_theo=D(129:192,6);
A_P_exp=D(193:256,6);
A_T_exp=D(257:320,6);
P_T_exp=D(321:384,6);
m=1;
for j=1:8
    for i=1:8
        A_P_mat_theo(i,j)=A_P_theo(m);
        A_T_mat_theo(i,j)=A_T_theo(m);
        P_T_mat_theo(i,j)=P_T_theo(m);
        A_P_mat_exp(i,j)=A_P_exp(m);
        A_T_mat_exp(i,j)=A_T_exp(m);
        P_T_mat_exp(i,j)=P_T_exp(m);
        m=m+1;
    end
end
%% djfjk

d=figure; 
subplot(3,2,1);
hm=heatmap((A_P_mat_theo));
caxis([-20 20])
hm.Colormap=turbo;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","1.22","4.88","19.53","78.125","312.5","1250","5000"];
% hm.Title="Abemaciclib (Y) vs PD0325901(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('Abem-PD-zip-theo.tiff','-dtiff',['-r' num2str(res)]);

% d=figure; 
subplot(3,2,2);
hm=heatmap((A_P_mat_exp));
caxis([-20 20])
hm.Colormap=turbo;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","1.22","4.88","19.53","78.125","312.5","1250","5000"];
% hm.Title="Abemaciclib (Y) vs PD0325901(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('Abem-PD-zip-exp.tiff','-dtiff',['-r' num2str(res)]);

subplot(3,2,3);
hm=heatmap((A_T_mat_theo));
caxis([-20 20])
hm.Colormap=turbo;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","0.0122","0.0488","0.1953","0.7813","3.125","12.50","50"];
% hm.Title="Abemaciclib (Y) vs TAK-960(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('Abem-TAK-zip-theo.tiff','-dtiff',['-r' num2str(res)]);

% d=figure;
subplot(3,2,4);
hm=heatmap((A_T_mat_exp));
caxis([-20 20])
hm.Colormap=turbo;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","0.0122","0.0488","0.1953","0.7813","3.125","12.50","50"];
% hm.Title="Abemaciclib (Y) vs TAK-960(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('Abem-TAK-zip-exp.tiff','-dtiff',['-r' num2str(res)]);

subplot(3,2,5);
hm=heatmap((P_T_mat_theo));
caxis([-20 20])
hm.Colormap=turbo;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","0.0122","0.0488","0.1953","0.7813","3.125","12.50","50"];
% hm.Title="PD0325901(Y) vs TAK-960(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('PD-TAK-zip-theo.tiff','-dtiff',['-r' num2str(res)]);

subplot(3,2,6);
hm=heatmap((P_T_mat_exp));
caxis([-20 20])
hm.Colormap=turbo;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","0.0122","0.0488","0.1953","0.7813","3.125","12.50","50"];
% hm.Title="PD0325901(Y) vs TAK-960(X)";
hm.FontSize=7;
res = 100;
set(gcf,'paperunits','inches','PaperPosition',[0 0 5 5]);
print('4D.tiff','-dtiff',['-r' num2str(res)]);
