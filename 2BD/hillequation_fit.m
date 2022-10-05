function Z = hillequation_fit(Fy,Fx,Dy,Dx)
Z=zeros(8,8);
for i = 1:8
    for j=1:8
        a = (Fy{i}(1)+Fy{i}(2)*(Dy(j)/Fy{i}(3))^Fy{i}(4))/(1+(Dy(j)/Fy{i}(3))^Fy{i}(4)); %Going downhill.
        b = (Fx{j}(1)+Fx{j}(2)*(Dx(i)/Fx{j}(3))^Fx{j}(4))/(1+(Dx(i)/Fx{j}(3))^Fx{j}(4)); %Going downhill.
        Zy = a;
        Zx = b;
        Z(j,i)=(Zx+Zy)/2;
    end
end
end