function Cell_response=Combination_total (params,steps,y0,Da,Dp,Dt,Mo)
% Mo=[params(7),params(8),params(9)];

% death parameters
MD=zeros(1,3);
Y = evalmarkovdata(steps,y0,Mo,MD);

%for Abemaciclib(M(2)) and PD0325901(M(1))
k=1;
for m=8:-1:1
    for n=1:8
        %inhibition
        MP12=Mo;
        MP12(2) = Mo(2)*(1-((Da(m)^params(4)))/((Da(m)^params(4))+(params(3)^params(4))));
        MP12(1) = Mo(1)*(1-((Dp(n)^params(2)))/((Dp(n)^params(2))+(params(1)^params(2))));
        
        %death
        MD=zeros(1,3);
        MD(1)=params(8)*Dp(n)/((Dp(n)+params(7)));
        MD(2)=params(10)*Da(m)/((Da(m)+params(9)));
        
        YP1P2{k,n}= evalmarkovdata(steps,y0,MP12,MD);
    end
    k=k+1;
end
%for Abemaciclib(M(2)) and TAK960(M(3))
k=1;
for m=8:-1:1
    for n=1:8
        %inhibition
        MP23=Mo;
        MP23(2) = Mo(2)*(1-((Da(m)^params(4)))/((Da(m)^params(4))+(params(3)^params(4))));
        MP23(3) = Mo(3)*(1-((Dt(n)^params(6)))/((Dt(n)^params(6))+(params(5)^params(6))));
        %death
        MD=zeros(1,3);
        MD(2)=params(10)*Da(m)/((Da(m)+params(9)));
        MD(3)=params(12)*Dt(n)/((Dt(n)+params(11)));
        
        YP2P3{k,n}= evalmarkovdata(steps,y0,MP23,MD);
    end
    k=k+1;
end

%for PD0325901(M(1)) and TAK960(M(3))
k=1;
for m=8:-1:1
    for n=1:8
        %inhibition
        MP13=Mo;
        MP13(1) = Mo(1)*(1-((Dp(m)^params(2)))/((Dp(m)^params(2))+(params(1)^params(2))));
        MP13(3) = Mo(3)*(1-((Dt(n)^params(6)))/((Dt(n)^params(6))+(params(5)^params(6))));
        %death
        MD=zeros(1,3);
        MD(1)=params(8)*Dp(m)/((Dp(m)+params(7)));
        MD(3)=params(12)*Dt(n)/((Dt(n)+params(11)));
        YP1P3{k,n}= evalmarkovdata(steps,y0,MP13,MD);
    end
    k=k+1;
end

function Y = evalmarkovdata(steps,y0,M,MD)
    T(1)=0;
    Y(1,1:3)=y0;
    for i=2:steps  
    Y(i,1)=Y(i-1,1)*(1-M(1)-MD(1))+Y(i-1,3)*M(3)*2;
    Y(i,2)=Y(i-1,1)*M(1)+Y(i-1,2)*(1-M(2)-MD(2));
    Y(i,3)=Y(i-1,2)*M(2)+Y(i-1,3)*(1-M(3)-MD(3));
    end
    T(i)=i-1;
end
for m=1:8
    for n=1:8
    Cell_response{1}(m,n)=sum(YP1P2{m,n}(steps,:))/sum(Y(steps,:));
    Cell_response{2}(m,n)=sum(YP2P3{m,n}(steps,:))/sum(Y(steps,:));
    Cell_response{3}(m,n)=sum(YP1P3{m,n}(steps,:))/sum(Y(steps,:));
    Cell_response{4}(m,n)=sum(YP1P2{m,n}(steps,:));
    Cell_response{5}(m,n)=sum(YP2P3{m,n}(steps,:));
    Cell_response{6}(m,n)=sum(YP1P3{m,n}(steps,:));

    % Cell_response(2)=sum(YP2(steps,:))/sum(Y(steps,:));
% Cell_response(3)=sum(YP3(steps,:))/sum(Y(steps,:));
    end
end
end
