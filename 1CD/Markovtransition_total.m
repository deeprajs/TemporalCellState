function [Cell_response,Y]=Markovtransition_total (params,steps,y0,D,Mo)

% death parameters
MD=zeros(1,3);
Y = evalmarkovdata(steps,y0,Mo,MD);

MP1=Mo;
MP1(1) = Mo(1)*(1-((D^params(2)))/((D^params(2))+(params(1)^params(2))));
MD=zeros(1,3);
MD(1)=params(8)*D/((D+params(7)));

YP1 = evalmarkovdata(steps,y0,MP1,MD);

MP2=Mo;
MP2(2) = Mo(2)*(1-((D^params(4)))/((D^params(4))+(params(3)^params(4))));
MD=zeros(1,3);
MD(2)=params(10)*D/((D+params(9)));
YP2 = evalmarkovdata(steps,y0,MP2,MD);


MP3=Mo;
MP3(3) = Mo(3)*(1-((D^params(6)))/((D^params(6))+(params(5)^params(6))));
MD=zeros(1,3);
MD(3)=params(12)*D/((D+params(11)));
YP3 = evalmarkovdata(steps,y0,MP3,MD);

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

Cell_response(1)=sum(YP1(steps,:))/sum(Y(steps,:));
Cell_response(2)=sum(YP2(steps,:))/sum(Y(steps,:));
Cell_response(3)=sum(YP3(steps,:))/sum(Y(steps,:));
Cell_response(4)=MP1(1);
Cell_response(5)=MP2(1);
Cell_response(6)=MP2(2);
Cell_response(7)=MP3(3);

end
