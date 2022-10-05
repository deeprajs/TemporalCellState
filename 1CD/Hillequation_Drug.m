function F = Hillequation_Drug(Ymax,EC50,N,D)

% The hill function used by lsqcurvefit.m.
% You get an n-value, Ymax, and IC50 from this function.
% x is a vector containing inital guess for the three x components.
% x(1) = n; x(2) = Ymax; x(3) = IC50

F = zeros(length(D),1);

for i = 1:length(D)
    a = (1-(Ymax*(D(i)^N))/((EC50^N)+(D(i)^N))); %Going downhill.
    F(i,1) = a;
end

end