%% 1C
% Get single drug dose response data for U87
[DR_data,~, ID]=xlsread("C:\Users\Deepraj\Documents\Dose Response\U87\GR_TagBFP+RFP_U87_significant_KIs_QC_r4.xlsx");
Abem_DR=DR_data(1:54,5)./DR_data(1:54,6);
PD_DR=DR_data(487:540,5)./DR_data(487:540,6);
TAK_DR=DR_data(649:702,5)./DR_data(649:702,6);

% Obtain initial fits for single drug dose responses
options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
Fa= lsqcurvefit(@hillequation_inh,[1;1;1;1],DR_data(1:54,4),1-Abem_DR,[0,0,0,0],[0,1,1,10],options);
Fp= lsqcurvefit(@hillequation_inh,[1;1;1;1],DR_data(1:54,4),1-PD_DR,[0,0,0,0],[0,1,1,10],options);
Ft= lsqcurvefit(@hillequation_inh,[1;1;1;1],DR_data(1:54,4),1-TAK_DR,[0,0,0,0],[0,1,1,10],options);

% Use initial fits to get single drug response parameters 
% for temporal cell state model

Ymax=[Fa(2),Fp(2),Ft(2)];
EC50=[Fa(3),Fp(3),Ft(3)];
N=[Fa(4),Fp(4),Ft(4)];
D=[0.001,0.00316,0.01,0.0316,0.1,0.316,1,3.16,10];
M=[0.05,0.14,0.14];
multistart=5;
for n=1:multistart
    [fmin(n,:),finSSE(n,1)]= estimatemarkovparams(Ymax,EC50,N,D,M);
n
end
[~,I] = min(finSSE);

D_fine=zeros(1,20);
D_fine(25)=10;
for i=24:-1:1
D_fine(i)=D_fine(i+1)/1.5;
end
steps=73;
y0=[64,19,17];
for i=1:length(D_fine)
    [GuessData(i,:),States]=Markovtransition_total (fmin(I,:),steps,y0,D_fine(i),M);
end


Abem_DR_e=std(reshape(Abem_DR,9,6),0,2)/2.45;
PD_DR_e=std(reshape(PD_DR,9,6),0,2)/2.45;
TAK_DR_e=std(reshape(TAK_DR,9,6),0,2)/2.45;

Abem_DR_m=mean(reshape(Abem_DR,9,6),2);
PD_DR_m=mean(reshape(PD_DR,9,6),2);
TAK_DR_m=mean(reshape(TAK_DR,9,6),2);

D=D*1000;
D_fine=D_fine*1000;
figure;
subplot(1,3,1);
hold on
set(gca, 'XScale', 'log');
% semilogx((D),P,'g.','MarkerSize', 15);
semilogx((D_fine),GuessData(:,1),'g-', 'MarkerSize', 15);
errorbar((D),PD_DR_m,PD_DR_e,'g.');
ylim([0 1.2])
xlim([1 10000])
title({"PD0325901", "MEK1/2 Inhibitor"},'FontSize',9)


subplot(1,3,2);
hold on
set(gca, 'XScale', 'log');
% semilogx((D),A,'b.','MarkerSize', 15);
semilogx((D_fine),GuessData(:,2),'b-', 'MarkerSize', 15);
errorbar((D),Abem_DR_m,Abem_DR_e,'b.');
ylim([0 1.2])
xlim([1 10000])
title({"Abemaciclib","CDK4/6 Inhibitor"},'FontSize',9)

subplot(1,3,3);
hold on
set(gca, 'XScale', 'log');
% semilogx((D),T,'r.','MarkerSize', 15);
semilogx((D_fine),GuessData(:,3),'r-', 'MarkerSize', 15);
errorbar((D),TAK_DR_m,TAK_DR_e,'r.');
ylim([0 1.2])
xlim([1 10000])
title({"TAK-960","PLK1 Inhibitor"},'FontSize',9)
res = 300;

set(gcf,'paperunits','inches','PaperPosition',[0 0 6 2]);
% set(gcf,'paperunits','inches','PaperPosition',[0 0 3.5 1.5]);
print('2B.tiff','-dtiff',['-r' num2str(res)]);

figure
hold on
yyaxis left
plot(0:1:72,sum(States,2),'k--');
% plot(0:1:72,States(:,1),'g--');
% plot(0:1:72,States(:,2),'b--');
% plot(0:1:72,States(:,3),'r--');
ylim([100 600])

yyaxis right
% plot(0:1:72,sum(States,2)./sum(States,2),'k-');
plot(0:1:72,States(:,1)./sum(States,2),'g-');
plot(0:1:72,States(:,2)./sum(States,2),'b-');
plot(0:1:72,States(:,3)./sum(States,2),'r-');
xlim([0 72])
ylim([0 1.1])
set(gcf,'paperunits','inches','PaperPosition',[0 0 2 2]);
print('2D.tiff','-dtiff',['-r' num2str(res)]);

%% 2C
clear;

[DR_data,~, ID]=xlsread("C:\Users\Deepraj\Documents\Dose Response\U251\GR_TagBFP+RFP_U251_All_KIs_QC_r3.csv");
Abem_DR=DR_data(1:54,5)./DR_data(1:54,6);
PD_DR=DR_data(343:396,5)./DR_data(343:396,6);
TAK_DR=DR_data(478:531,5)./DR_data(478:531,6);

Abem_DR_e=std(reshape(Abem_DR,9,6),0,2)/2.45;
PD_DR_e=std(reshape(PD_DR,9,6),0,2)/2.45;
TAK_DR_e=std(reshape(TAK_DR,9,6),0,2)/2.45;

Abem_DR_m=mean(reshape(Abem_DR,9,6),2);
PD_DR_m=mean(reshape(PD_DR,9,6),2);
TAK_DR_m=mean(reshape(TAK_DR,9,6),2);

options = optimoptions('lsqcurvefit','TolFun',1e-14,'TolX',1e-14);
Fa= lsqcurvefit(@hillequation_inh,[1;1;1;1],DR_data(1:54,4),1-Abem_DR,[0,0,0,0],[0,1,1,10],options);
Fp= lsqcurvefit(@hillequation_inh,[1;1;1;1],DR_data(1:54,4),1-PD_DR,[0,0,0,0],[0,1,1,10],options);
Ft= lsqcurvefit(@hillequation_inh,[1;1;1;1],DR_data(1:54,4),1-TAK_DR,[0,0,0,0],[0,1,1,10],options);

Ymax=[Fa(2),Fp(2),Ft(2)];
EC50=[Fa(3),Fp(3),Ft(3)];
N=[Fa(4),Fp(4),Ft(4)];
D=[0.001,0.00316,0.01,0.0316,0.1,0.316,1,3.16,10];
M=[0.075,0.16,0.16]; 
% [fmin,finSSE]= estimatemarkovparams(Ymax,EC50,N,D);

multistart=5;
for n=1:multistart
    [fmin(n,:),finSSE(n,1)]= estimatemarkovparams(Ymax,EC50,N,D,M);
n
end
D_fine=zeros(1,20);
D_fine(25)=10;
for i=24:-1:1
D_fine(i)=D_fine(i+1)/1.5;
end
steps=73;
y0=[58,23,19];
% fmin(1)=0.07;
for i=1:length(D_fine)
    [GuessData(i,:),States]=Markovtransition_total (fmin(1,:),steps,y0,D_fine(i),M);
end
D=D*1000;
D_fine=D_fine*1000;
figure;
subplot(1,3,1);
hold on
set(gca, 'XScale', 'log');
% semilogx((D),P,'g.','MarkerSize', 15);
semilogx((D_fine),GuessData(:,1),'g-', 'MarkerSize', 15);
errorbar((D),PD_DR_m,PD_DR_e,'g.');
ylim([0 1.2])
xlim([1 10000])
title({"PD0325901", "MEK1/2 Inhibitor"},'FontSize',9)

subplot(1,3,2);
hold on
set(gca, 'XScale', 'log');
% semilogx((D),A,'b.','MarkerSize', 15);
semilogx((D_fine),GuessData(:,2),'b-', 'MarkerSize', 15);
errorbar((D),Abem_DR_m,Abem_DR_e,'b.');
ylim([0 1.2])
xlim([1 10000])
title({"Abemaciclib","CDK4/6 Inhibitor"},'FontSize',9)

subplot(1,3,3);
hold on
set(gca, 'XScale', 'log');
% semilogx((D),T,'r.','MarkerSize', 15);
semilogx((D_fine),GuessData(:,3),'r-', 'MarkerSize', 15);
errorbar((D),TAK_DR_m,TAK_DR_e,'r.');
ylim([0 1.2])
xlim([1 10000])
title({"TAK-960","PLK1 Inhibitor"},'FontSize',9)
res = 300;

set(gcf,'paperunits','inches','PaperPosition',[0 0 6 2]);
% set(gcf,'paperunits','inches','PaperPosition',[0 0 3.5 1.5]);
print('2C.tiff','-dtiff',['-r' num2str(res)]);
figure
hold on
yyaxis left
plot(0:1:72,sum(States,2),'k--');
% plot(0:1:72,States(:,1),'g--');
% plot(0:1:72,States(:,2),'b--');
% plot(0:1:72,States(:,3),'r--');
ylim([100 1000])

yyaxis right
% plot(0:1:72,sum(States,2)./sum(States,2),'k-');
plot(0:1:72,States(:,1)./sum(States,2),'g-');
plot(0:1:72,States(:,2)./sum(States,2),'b-');
plot(0:1:72,States(:,3)./sum(States,2),'r-');
xlim([0 72])
ylim([0 1.1])
set(gcf,'paperunits','inches','PaperPosition',[0 0 2 2]);
print('2E.tiff','-dtiff',['-r' num2str(res)]);