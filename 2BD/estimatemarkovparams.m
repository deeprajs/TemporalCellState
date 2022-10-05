function [fmin,finSSE]= estimatemarkovparams(Ymax_D,EC50_D,N_D,D,M)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

conoptions = optimoptions('fmincon','MaxFunctionEvaluations',10000);

A=Hillequation_Drug(Ymax_D(1),EC50_D(1),N_D(1),D);
P=Hillequation_Drug(Ymax_D(2),EC50_D(2),N_D(2),D);
T=Hillequation_Drug(Ymax_D(3),EC50_D(3),N_D(3),D);

% DoseResponse(1:length(D),1)=A;
% DoseResponse(1:length(D),2)=P;
% DoseResponse(1:length(D),3)=T;

% Set upper and lower bounds for edge weights
%          lb = [0 0 0 0 0 0 0.02 0.1 0.1];
%          ub = [1 10 1 10 1 10 0.08 0.2 0.2];
         lb = [0 0 0 0 0 0 0 0 0 0 0 0];
         ub = [1 10 1 10 1 10 1000 1 1000 1 1000 0.02];
         
initial=[rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1)];
% initial=[rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1),rand(1)];

% initial=[EC50_D(1),N_D(1),EC50_D(2),N_D(2),EC50_D(3),N_D(3)];
[fmin,finSSE]=fmincon(@nestedfun,initial,[],[],[],[],lb,ub,[],conoptions);
         
        
function SSE = nestedfun(r2)
        
    % Solve the ODEs 
steps=73;
y0=[58,23,19];

for i=1:length(D)
    Guess(i,:)=Markovtransition_total (r2,steps,y0,D(i),M);
end

%  SSE = (sum(sum(((A-Guess(:,1))./A).^2)+sum(((P-Guess(:,2))./P).^2)+sum(((T-Guess(:,3))./T).^2))).^0.5;
 SSE = (sum(sum(((P-Guess(:,1))).^2)+sum(((A-Guess(:,2))).^2)+sum(((T-Guess(:,3))).^2))).^0.5;
         
end
end

