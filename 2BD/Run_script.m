Ymax=[0.952,1,0.959];
EC50=[0.096,0.34,0.003];
N=[0.781,0.551,1.653];
D=[0.001,0.00316,0.01,0.0316,0.1,0.316,1,3.16,10];
M=[0.075,0.16,0.16]; 
% [fmin,finSSE]= estimatemarkovparams(Ymax,EC50,N,D);

multistart=1;
for n=1:multistart
    [fmin(n,:),finSSE(n,1)]= estimatemarkovparams(Ymax,EC50,N,D,M);
n
end

%% Combination data
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
Mo=[0.075,0.16,0.16]; 
steps=73;
y0=[58,23,19];
Combinationdata=Combination_total (fmin(1,:),steps,y0,Da,Dp,Dt,Mo);


%% Abem-PD theo Dose response
Counts_theo=Combinationdata{1,1}/Combinationdata{1,1}(8,1);
d=figure; 
subplot(3,2,1);
hm=heatmap(Counts_theo);
caxis([0 1])
hm.Colormap=hot;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","1.22","4.88","19.53","78.125","312.5","1250","5000"];
% hm.Title="Abemaciclib (Y) vs PD0325901(X)";
% % hm.Title="Drug-1 (Y) vs Drug-2 (X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('Abem-PD-theo.tiff','-dtiff',['-r' num2str(res)]);
% %% Abem-PD Exp dose response
A=xlsread('E:\U251\deepraj_u251_abem_PD_bio1\MyExpt_Image.csv');
B=xlsread('E:\U251\deepraj_u251_abem_PD_bio2\MyExpt_Image.csv');
C=xlsread('E:\U251\deepraj_u251_abem_PD_bio3\MyExpt_Image.csv');

for i=1:64
    Counts_1(i,1)=A(1+8*(i-1))+A(2+8*(i-1))+A(3+8*(i-1))+A(4+8*(i-1))-A(5+8*(i-1))-A(6+8*(i-1))-A(7+8*(i-1))-A(8+8*(i-1));
    Counts_2(i,1)=B(1+8*(i-1))+B(2+8*(i-1))+B(3+8*(i-1))+B(4+8*(i-1))-B(5+8*(i-1))-B(6+8*(i-1))-B(7+8*(i-1))-B(8+8*(i-1));
    Counts_3(i,1)=C(1+8*(i-1))+C(2+8*(i-1))+C(3+8*(i-1))+C(4+8*(i-1))-C(5+8*(i-1))-C(6+8*(i-1))-C(7+8*(i-1))-C(8+8*(i-1));

end
Counts_1_mat_AP=(reshape(Counts_1,[8,8]).');
Counts_1_mat_AP=Counts_1_mat_AP/Counts_1_mat_AP(8,1);
Counts_2_mat_AP=reshape(Counts_2,[8,8]).';
Counts_2_mat_AP=Counts_2_mat_AP/Counts_2_mat_AP(8,1);
Counts_3_mat_AP=reshape(Counts_3,[8,8]).';
Counts_3_mat_AP=Counts_3_mat_AP/Counts_3_mat_AP(8,1);
Counts_tot=Counts_1_mat_AP+Counts_2_mat_AP+Counts_3_mat_AP;
% GC=Counts_tot-37.5*3*5;
Counts_tot_AP=Counts_tot/Counts_tot(8,1);
% d=figure; 
subplot(3,2,2);
hm=heatmap(Counts_tot_AP);
caxis([0 1])
hm.Colormap=hot;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","1.22","4.88","19.53","78.125","312.5","1250","5000"];
% hm.Title="Abemaciclib (Y) vs PD0325901(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('Abem-PD-exp.tiff','-dtiff',['-r' num2str(res)]);
% save('Abem_PD_exp.mat','Counts_tot')

% %% Abem-Tak
Counts_theo=Combinationdata{1,2}/Combinationdata{1,2}(8,1);
% d=figure; 
subplot(3,2,3);
hm=heatmap(Counts_theo);
caxis([0 1])
hm.Colormap=hot;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","0.0122","0.0488","0.1953","0.7813","3.125","12.50","50"];
% hm.Title="Abemaciclib (Y) vs TAK-960(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('Abem-TAK-theo.tiff','-dtiff',['-r' num2str(res)]);

% %% Abem TAK exp
A=xlsread('E:\U251\deepraj_abem_TAK_U251_bio1\MyExpt_Image.csv');
B=xlsread('E:\U251\deepraj_abem_TAK_U251_bio2\MyExpt_Image.csv');
C=xlsread('E:\U251\deepraj_abem_TAK_U251_bio3\MyExpt_Image.csv');

for i=1:64
    Counts_1(i,1)=A(1+8*(i-1))+A(2+8*(i-1))+A(3+8*(i-1))+A(4+8*(i-1))-A(5+8*(i-1))-A(6+8*(i-1))-A(7+8*(i-1))-A(8+8*(i-1));
    Counts_2(i,1)=B(1+8*(i-1))+B(2+8*(i-1))+B(3+8*(i-1))+B(4+8*(i-1))-B(5+8*(i-1))-B(6+8*(i-1))-B(7+8*(i-1))-B(8+8*(i-1));
    Counts_3(i,1)=C(1+8*(i-1))+C(2+8*(i-1))+C(3+8*(i-1))+C(4+8*(i-1))-C(5+8*(i-1))-C(6+8*(i-1))-C(7+8*(i-1))-C(8+8*(i-1));

end
Counts_1_mat_AT=reshape(Counts_1,[8,8]).';
Counts_1_mat_AT=Counts_1_mat_AT/Counts_1_mat_AT(8,1);
Counts_2_mat_AT=reshape(Counts_2,[8,8]).';
Counts_2_mat_AT=Counts_2_mat_AT/Counts_2_mat_AT(8,1);
Counts_3_mat_AT=reshape(Counts_3,[8,8]).';
Counts_3_mat_AT=Counts_3_mat_AT/Counts_3_mat_AT(8,1);
Counts_tot=Counts_1_mat_AT+Counts_2_mat_AT+Counts_3_mat_AT;
% GC=Counts_tot-37.5*3*5;

Counts_tot_AT=Counts_tot/Counts_tot(8,1);
% d=figure; 
subplot(3,2,4);
hm=heatmap(Counts_tot_AT);
caxis([0 1])
hm.Colormap=hot;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","0.0122","0.0488","0.1953","0.7813","3.125","12.50","50"];
% hm.Title="Abemaciclib (Y) vs TAK-960(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('Abem-TAK-exp.tiff','-dtiff',['-r' num2str(res)]);
% save('Abem_TAK_exp.mat','Counts_tot')

% %% PD-Tak
Counts_theo=Combinationdata{1,3}/Combinationdata{1,3}(8,1);
% d=figure; 
subplot(3,2,5);
hm=heatmap(Counts_theo);
caxis([0 1])
hm.Colormap=hot;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","0.0122","0.0488","0.1953","0.7813","3.125","12.50","50"];
% hm.Title="PD0325901 (Y) vs TAK-960(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('PD-TAK-theo.tiff','-dtiff',['-r' num2str(res)]);

% %% PD-tak exp
A=xlsread('E:\U251\deepraj_PD_TAK_U251_bio1\MyExpt_Image.csv');
B=xlsread('E:\U251\deepraj_PD_TAK_U251_bio2\MyExpt_Image.csv');
C=xlsread('E:\U251\deepraj_PD_TAK_U251_bio3\MyExpt_Image.csv');

for i=1:64
    Counts_1(i,1)=A(1+8*(i-1))+A(2+8*(i-1))+A(3+8*(i-1))+A(4+8*(i-1))-A(5+8*(i-1))-A(6+8*(i-1))-A(7+8*(i-1))-A(8+8*(i-1));
    Counts_2(i,1)=B(1+8*(i-1))+B(2+8*(i-1))+B(3+8*(i-1))+B(4+8*(i-1))-B(5+8*(i-1))-B(6+8*(i-1))-B(7+8*(i-1))-B(8+8*(i-1));
    Counts_3(i,1)=C(1+8*(i-1))+C(2+8*(i-1))+C(3+8*(i-1))+C(4+8*(i-1))-C(5+8*(i-1))-C(6+8*(i-1))-C(7+8*(i-1))-C(8+8*(i-1));

end
Counts_1_mat_PT=reshape(Counts_1,[8,8]).';
Counts_1_mat_PT=Counts_1_mat_PT/Counts_1_mat_PT(8,1);
Counts_2_mat_PT=reshape(Counts_2,[8,8]).';
Counts_2_mat_PT=Counts_2_mat_PT/Counts_2_mat_PT(8,1);
Counts_3_mat_PT=reshape(Counts_3,[8,8]).';
Counts_3_mat_PT=Counts_3_mat_PT/Counts_3_mat_PT(8,1);
% GC=Counts_tot-37.5*3*5;

Counts_tot_PT=Counts_tot/Counts_tot(8,1);
% d=figure; 
subplot(3,2,6);
hm=heatmap(Counts_tot_PT);
caxis([0 1])
hm.Colormap=hot;
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% hm.YData = flip(["0","1.22","4.88","19.53","78.125","312.5","1250","5000"]);
% hm.XData = ["0","0.0122","0.0488","0.1953","0.7813","3.125","12.50","50"];
% hm.Title="PD0325901 (Y) vs TAK-960(X)";
hm.FontSize=7;
res = 100;
% set(gcf,'paperunits','inches','PaperPosition',[0 0 2 1.5]);
% print('PD-TAK-exp.tiff','-dtiff',['-r' num2str(res)]);
% save('PD_TAK_exp.mat','Counts_tot')
set(gcf,'paperunits','inches','PaperPosition',[0 0 4.5 4.5]);
print('3B.tiff','-dtiff',['-r' num2str(res)]);
%% EOB
m=1;
Inh=(1-Combinationdata{1,m});
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dp,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Fp,flip(Da),Dp);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AP_theo=Inh_fit-Bliss;

    m=1;
Inh=1-Counts_tot_AP;
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dp,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Fp,flip(Da),Dp);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AP_exp=Inh_fit-Bliss;
    
    EOB_AP(1,1)=mean(mean(EOB_AP_theo(1:3,6:8)));
    EOB_AP(2,1)=mean(mean(EOB_AP_theo(1:3,2:5)));
    EOB_AP(3,1)=mean(mean(EOB_AP_theo(4:7,6:8)));
    EOB_AP(4,1)=mean(mean(EOB_AP_theo(4:7,2:5)));
    
    EOB_AP(1,2)=mean(mean(EOB_AP_exp(1:3,6:8)));
    EOB_AP(2,2)=mean(mean(EOB_AP_exp(1:3,2:5)));
    EOB_AP(3,2)=mean(mean(EOB_AP_exp(4:7,6:8)));
    EOB_AP(4,2)=mean(mean(EOB_AP_exp(4:7,2:5)));
    
%     subplot(3,1,1);
%     bar(EOB*100);
%     ylim([-20 10])
%         set(gca,'xticklabel',{[]})

    %AT
m=2;
Inh=(1-Combinationdata{1,m});
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Ft,flip(Da),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AT_theo=Inh_fit-Bliss;

    m=1;
Inh=1-Counts_tot_AT;
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Ft,flip(Da),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AT_exp=Inh_fit-Bliss;
    
    EOB_AT(1,1)=mean(mean(EOB_AT_theo(1:3,6:8)));
    EOB_AT(2,1)=mean(mean(EOB_AT_theo(1:3,2:5)));
    EOB_AT(3,1)=mean(mean(EOB_AT_theo(4:7,6:8)));
    EOB_AT(4,1)=mean(mean(EOB_AT_theo(4:7,2:5)));
    
    EOB_AT(1,2)=mean(mean(EOB_AT_exp(1:3,6:8)));
    EOB_AT(2,2)=mean(mean(EOB_AT_exp(1:3,2:5)));
    EOB_AT(3,2)=mean(mean(EOB_AT_exp(4:7,6:8)));
    EOB_AT(4,2)=mean(mean(EOB_AT_exp(4:7,2:5)));
    
%     subplot(3,1,2);
%     bar(EOB*100);
%     ylim([-20 10])
%         set(gca,'xticklabel',{[]})

       %PT
m=3;
Inh=(1-Combinationdata{1,m});
%     Inh=1-Counts_tot_AP;
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fp,Ft,flip(Da),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_PT_theo=Inh_fit-Bliss;

    m=1;
Inh=1-Counts_tot_PT;
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Ft,flip(Da),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_PT_exp=Inh_fit-Bliss;
    
    EOB_PT(1,1)=mean(mean(EOB_PT_theo(1:3,6:8)));
    EOB_PT(2,1)=mean(mean(EOB_PT_theo(1:3,2:5)));
    EOB_PT(3,1)=mean(mean(EOB_PT_theo(4:7,6:8)));
    EOB_PT(4,1)=mean(mean(EOB_PT_theo(4:7,2:5)));
    
    EOB_PT(1,2)=mean(mean(EOB_PT_exp(1:3,6:8)));
    EOB_PT(2,2)=mean(mean(EOB_PT_exp(1:3,2:5)));
    EOB_PT(3,2)=mean(mean(EOB_PT_exp(4:7,6:8)));
    EOB_PT(4,2)=mean(mean(EOB_PT_exp(4:7,2:5)));
    
    
    %% erorbars
Inh=1-Counts_1_mat_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dp,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Fp,flip(Da),Dp);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AP_exp=Inh_fit-Bliss;
    
    EOB_AP_err(1,1)=mean(mean(EOB_AP_exp(1:3,6:8)));
    EOB_AP_err(2,1)=mean(mean(EOB_AP_exp(1:3,2:5)));
    EOB_AP_err(3,1)=mean(mean(EOB_AP_exp(4:7,6:8)));
    EOB_AP_err(4,1)=mean(mean(EOB_AP_exp(4:7,2:5)));
 
%     Counts_2_mat_AP(8,:)=Counts_tot_AP(8,:);
%     Counts_2_mat_AP(:,1)=Counts_tot_AP(:,1);

Inh=1-Counts_2_mat_AP;
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dp,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Fp,flip(Da),Dp);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AP_exp=Inh_fit-Bliss;
    
    EOB_AP_err(1,2)=mean(mean(EOB_AP_exp(1:3,6:8)));
    EOB_AP_err(2,2)=mean(mean(EOB_AP_exp(1:3,2:5)));
    EOB_AP_err(3,2)=mean(mean(EOB_AP_exp(4:7,6:8)));
    EOB_AP_err(4,2)=mean(mean(EOB_AP_exp(4:7,2:5)));
    
    Inh=1-Counts_3_mat_AP;
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dp,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Fp,flip(Da),Dp);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AP_exp=Inh_fit-Bliss;
    
    EOB_AP_err(1,3)=mean(mean(EOB_AP_exp(1:3,6:8)));
    EOB_AP_err(2,3)=mean(mean(EOB_AP_exp(1:3,2:5)));
    EOB_AP_err(3,3)=mean(mean(EOB_AP_exp(4:7,6:8)));
    EOB_AP_err(4,3)=mean(mean(EOB_AP_exp(4:7,2:5)));
    
    Inh=1-Counts_1_mat_AT;
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Ft,flip(Da),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AT_exp=Inh_fit-Bliss;

    EOB_AT_err(1,1)=mean(mean(EOB_AT_exp(1:3,6:8)));
    EOB_AT_err(2,1)=mean(mean(EOB_AT_exp(1:3,2:5)));
    EOB_AT_err(3,1)=mean(mean(EOB_AT_exp(4:7,6:8)));
    EOB_AT_err(4,1)=mean(mean(EOB_AT_exp(4:7,2:5)));
    
        Inh=1-Counts_2_mat_AT;
%     Inh=1-Counts_tot_AP;
Da=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Ft,flip(Dp),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AT_exp=Inh_fit-Bliss;

    EOB_AT_err(1,2)=mean(mean(EOB_AT_exp(1:3,6:8)));
    EOB_AT_err(2,2)=mean(mean(EOB_AT_exp(1:3,2:5)));
    EOB_AT_err(3,2)=mean(mean(EOB_AT_exp(4:7,6:8)));
    EOB_AT_err(4,2)=mean(mean(EOB_AT_exp(4:7,2:5)));
    
    Inh=1-Counts_3_mat_AT;
%     Inh=1-Counts_tot_AP;
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fa{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fa,Ft,flip(Dp),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_AT_exp=Inh_fit-Bliss;

    EOB_AT_err(1,3)=mean(mean(EOB_AT_exp(1:3,6:8)));
    EOB_AT_err(2,3)=mean(mean(EOB_AT_exp(1:3,2:5)));
    EOB_AT_err(3,3)=mean(mean(EOB_AT_exp(4:7,6:8)));
    EOB_AT_err(4,3)=mean(mean(EOB_AT_exp(4:7,2:5)));
    
    Inh=1-Counts_1_mat_PT;
%     Inh=1-Counts_tot_AP;
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fp,Ft,flip(Dp),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_PT_exp=Inh_fit-Bliss;

    EOB_PT_err(1,1)=mean(mean(EOB_PT_exp(1:3,6:8)));
    EOB_PT_err(2,1)=mean(mean(EOB_PT_exp(1:3,2:5)));
    EOB_PT_err(3,1)=mean(mean(EOB_PT_exp(4:7,6:8)));
    EOB_PT_err(4,1)=mean(mean(EOB_PT_exp(4:7,2:5)));
    
        Inh=1-Counts_2_mat_PT;
%     Inh=1-Counts_tot_AP;
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fp,Ft,flip(Dp),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_PT_exp=Inh_fit-Bliss;

    EOB_PT_err(1,2)=mean(mean(EOB_PT_exp(1:3,6:8)));
    EOB_PT_err(2,2)=mean(mean(EOB_PT_exp(1:3,2:5)));
    EOB_PT_err(3,2)=mean(mean(EOB_PT_exp(4:7,6:8)));
    EOB_PT_err(4,2)=mean(mean(EOB_PT_exp(4:7,2:5)));
    
    Inh=1-Counts_3_mat_PT;
%     Inh=1-Counts_tot_AP;
Dp=[0,0.00122,0.00488,0.01953,0.078125,0.3125,1.250,5];
Dt=[0,0.0000122,0.0000488,0.0001953,0.0007813,0.003125,0.0125,0.050];
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
for i=1:8
Fp{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Da,flip(Inh(:,i)),[0,0,0,0],[Inh(8,i),1,1,10],options);
Ft{i} = lsqcurvefit(@hillequation_inh,[1;1;1;1],Dt,(Inh(i,:))',[0,0,0,0],[Inh(1,i),1,1,10],options);
end
Inh_fit = hillequation_fit(Fp,Ft,flip(Dp),Dt);
 
    Bliss(8,:)=Inh_fit(8,:);
    Bliss(:,1)=Inh_fit(:,1);
    for i=1:7
        for j=2:8
        Bliss(i,j)=Inh_fit(i,1)+Inh_fit(8,j)-Inh_fit(i,1)*Inh_fit(8,j);
        end
    end
    EOB_PT_exp=Inh_fit-Bliss;

    EOB_PT_err(1,3)=mean(mean(EOB_PT_exp(1:3,6:8)));
    EOB_PT_err(2,3)=mean(mean(EOB_PT_exp(1:3,2:5)));
    EOB_PT_err(3,3)=mean(mean(EOB_PT_exp(4:7,6:8)));
    EOB_PT_err(4,3)=mean(mean(EOB_PT_exp(4:7,2:5)));
    %% dkls
figure
hold on   
area([-0.15,0.05],[0.05,0.05],'FaceColor','r');
alpha(0.1)
area([0.05,-0.15],[-0.15,-0.15],'FaceColor','b');
alpha(0.1)

plot(-0.15:0.2:0.05,-0.15:0.2:0.05,'k--');
plot(-0.15:0.2:0.05,zeros(2),'k--');
plot(zeros(2),-0.15:0.2:0.05,'k--');

fmt = {'p','o','^','s'};
col = {'g','b','r'};

f = 0;
hold on
for n = 1:4
scatter(EOB_AP(n,1),mean(EOB_AP_err(n,:),2),60,'g',[fmt{f+n}],'filled');
errorbar(EOB_AP(n,1),mean(EOB_AP_err(n,:),2),std(EOB_AP_err(n,:)/1.73,0,2),'g');
scatter(EOB_AT(n,1),mean(EOB_AT_err(n,:),2),60,'b',[fmt{f+n}],'filled');
errorbar(EOB_AT(n,1),mean(EOB_AT_err(n,:),2),std(EOB_AT_err(n,:)/1.73,0,2),'b');
scatter(EOB_PT(n,1),mean(EOB_PT_err(n,:),2),60,'r',[fmt{f+n}],'filled');
errorbar(EOB_PT(n,1),mean(EOB_PT_err(n,:),2),std(EOB_PT_err(n,:)/1.73,0,2),'r');

end
ylim([-0.15 0.05])
xlim([-0.15 0.05])
xlabel('Model','FontWeight','bold') 
ylabel('Experiment','FontWeight','bold') 
set(gcf,'paperunits','inches','PaperPosition',[0 0 4 3.6]);
print('3D.tiff','-dtiff',['-r' num2str(100)]);