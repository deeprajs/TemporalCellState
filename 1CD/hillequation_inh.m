function F = hillequation_inh(x,D)
% The hill function used by lsqcurvefit.m.
% You get an n-value, Ymax, and IC50 from this function.
% x is a vector containing inital guess for the three x components.
% x(1) = n; x(2) = Ymax; x(3) = IC50

F = zeros(length(D),1);
for i = 1:length(D)
    a = (x(1)+x(2)*(D(i)/x(3))^x(4))/(1+(D(i)/x(3))^x(4)); %Going downhill.
    F(i) = a;
end

end