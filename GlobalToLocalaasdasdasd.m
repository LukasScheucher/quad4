function xi=GlobalToLocal(elenodes,Xg)

x1=elenodes(1,1);
x2=elenodes(2,1);
x3=elenodes(3,1);
x4=elenodes(4,1);

y1=elenodes(1,2);
y2=elenodes(2,2);
y3=elenodes(3,2);
y4=elenodes(4,2);

X=Xg(1);

Y=Xg(2);


 %xi(1)=-(x1 - 4*X + x2 + x3 + x4 + 4*Y*x1 - 4*Y*x2 + 4*Y*x3 - 4*Y*x4 + (16*X^2 - 32*X*Y*x1 + 32*X*Y*x2 - 32*X*Y*x3 + 32*X*Y*x4 - 16*X*x1*y1^2 + 32*X*x1*y1*y3 + 32*X*x1*y1 + 16*X*x1*y2^2 - 32*X*x1*y2*y4 - 16*X*x1*y3^2 + 16*X*x1*y4^2 - 8*X*x1 + 16*X*x2*y1^2 - 32*X*x2*y1*y3 - 16*X*x2*y1 - 16*X*x2*y2^2 + 32*X*x2*y2*y4 - 16*X*x2*y2 + 16*X*x2*y3^2 - 16*X*x2*y3 - 16*X*x2*y4^2 + 16*X*x2*y4 - 8*X*x2 - 16*X*x3*y1^2 + 32*X*x3*y1*y3 + 16*X*x3*y2^2 - 32*X*x3*y2*y4 - 16*X*x3*y3^2 + 32*X*x3*y3 + 16*X*x3*y4^2 - 8*X*x3 + 16*X*x4*y1^2 - 32*X*x4*y1*y3 - 16*X*x4*y1 - 16*X*x4*y2^2 + 32*X*x4*y2*y4 + 16*X*x4*y2 + 16*X*x4*y3^2 - 16*X*x4*y3 - 16*X*x4*y4^2 - 16*X*x4*y4 - 8*X*x4 + 16*Y^2*x1^2 - 32*Y^2*x1*x2 + 32*Y^2*x1*x3 - 32*Y^2*x1*x4 + 16*Y^2*x2^2 - 32*Y^2*x2*x3 + 32*Y^2*x2*x4 + 16*Y^2*x3^2 - 32*Y^2*x3*x4 + 16*Y^2*x4^2 - 32*Y*x1^2*y3 - 8*Y*x1^2 + 16*Y*x1*x2*y1 - 16*Y*x1*x2*y2 + 48*Y*x1*x2*y3 + 16*Y*x1*x2*y4 - 32*Y*x1*x3*y1 - 32*Y*x1*x3*y3 + 48*Y*x1*x3 + 16*Y*x1*x4*y1 + 16*Y*x1*x4*y2 + 48*Y*x1*x4*y3 - 16*Y*x1*x4*y4 - 16*Y*x2^2*y1 + 16*Y*x2^2*y2 - 16*Y*x2^2*y3 - 16*Y*x2^2*y4 + 8*Y*x2^2 + 48*Y*x2*x3*y1 - 16*Y*x2*x3*y2 + 16*Y*x2*x3*y3 + 16*Y*x2*x3*y4 - 32*Y*x2*x4*y1 - 32*Y*x2*x4*y3 - 48*Y*x2*x4 - 32*Y*x3^2*y1 - 8*Y*x3^2 + 48*Y*x3*x4*y1 + 16*Y*x3*x4*y2 + 16*Y*x3*x4*y3 - 16*Y*x3*x4*y4 - 16*Y*x4^2*y1 - 16*Y*x4^2*y2 - 16*Y*x4^2*y3 + 16*Y*x4^2*y4 + 8*Y*x4^2 + 16*x1^2*y3^2 + 8*x1^2*y3 + x1^2 - 16*x1*x2*y1*y3 - 4*x1*x2*y1 + 16*x1*x2*y2*y3 + 4*x1*x2*y2 - 16*x1*x2*y3^2 - 16*x1*x2*y3*y4 + 4*x1*x2*y3 - 4*x1*x2*y4 + 2*x1*x2 + 16*x1*x3*y1^2 - 24*x1*x3*y1 - 16*x1*x3*y2^2 + 32*x1*x3*y2*y4 + 16*x1*x3*y3^2 - 24*x1*x3*y3 - 16*x1*x3*y4^2 + 2*x1*x3 - 16*x1*x4*y1*y3 - 4*x1*x4*y1 - 16*x1*x4*y2*y3 - 4*x1*x4*y2 - 16*x1*x4*y3^2 + 16*x1*x4*y3*y4 + 4*x1*x4*y3 + 4*x1*x4*y4 + 2*x1*x4 + 4*x2^2*y1^2 - 8*x2^2*y1*y2 + 8*x2^2*y1*y3 + 8*x2^2*y1*y4 - 4*x2^2*y1 + 4*x2^2*y2^2 - 8*x2^2*y2*y3 - 8*x2^2*y2*y4 + 4*x2^2*y2 + 4*x2^2*y3^2 + 8*x2^2*y3*y4 - 4*x2^2*y3 + 4*x2^2*y4^2 - 4*x2^2*y4 + x2^2 - 16*x2*x3*y1^2 + 16*x2*x3*y1*y2 - 16*x2*x3*y1*y3 - 16*x2*x3*y1*y4 + 4*x2*x3*y1 + 4*x2*x3*y2 - 4*x2*x3*y3 - 4*x2*x3*y4 + 2*x2*x3 - 8*x2*x4*y1^2 + 48*x2*x4*y1*y3 + 24*x2*x4*y1 + 8*x2*x4*y2^2 - 16*x2*x4*y2*y4 - 8*x2*x4*y3^2 + 24*x2*x4*y3 + 8*x2*x4*y4^2 + 2*x2*x4 + 16*x3^2*y1^2 + 8*x3^2*y1 + x3^2 - 16*x3*x4*y1^2 - 16*x3*x4*y1*y2 - 16*x3*x4*y1*y3 + 16*x3*x4*y1*y4 + 4*x3*x4*y1 - 4*x3*x4*y2 - 4*x3*x4*y3 + 4*x3*x4*y4 + 2*x3*x4 + 4*x4^2*y1^2 + 8*x4^2*y1*y2 + 8*x4^2*y1*y3 - 8*x4^2*y1*y4 - 4*x4^2*y1 + 4*x4^2*y2^2 + 8*x4^2*y2*y3 - 8*x4^2*y2*y4 - 4*x4^2*y2 + 4*x4^2*y3^2 - 8*x4^2*y3*y4 - 4*x4^2*y3 + 4*x4^2*y4^2 + 4*x4^2*y4 + x4^2)^(1/2) - 2*x1*y1 + 2*x1*y2 - 2*x1*y3 - 2*x3*y1 - 2*x1*y4 + 4*x2*y3 - 2*x3*y2 + 4*x4*y1 - 2*x3*y3 + 2*x3*y4)/(2*x2 + x1*(2*y1 - 2) - 2*x2*y1 + y4*(2*x1 - 2*x2 + 2*x3 - 2*x4) + x3*(2*y1 - 2*y2 - 2*y3 + 2) - x4*(2*y1 - 2*y2 - 2*y3 + 2) - y2*(2*x1 - 2*x2) - y3*(2*x1 - 2*x2));
 xi(1)=(4*X - x1 - x2 - x3 - x4 - 4*Y*x1 + 4*Y*x2 - 4*Y*x3 + 4*Y*x4 + (16*X^2 - 32*X*Y*x1 + 32*X*Y*x2 - 32*X*Y*x3 + 32*X*Y*x4 - 16*X*x1*y1^2 + 32*X*x1*y1*y3 + 32*X*x1*y1 + 16*X*x1*y2^2 - 32*X*x1*y2*y4 - 16*X*x1*y3^2 + 16*X*x1*y4^2 - 8*X*x1 + 16*X*x2*y1^2 - 32*X*x2*y1*y3 - 16*X*x2*y1 - 16*X*x2*y2^2 + 32*X*x2*y2*y4 - 16*X*x2*y2 + 16*X*x2*y3^2 - 16*X*x2*y3 - 16*X*x2*y4^2 + 16*X*x2*y4 - 8*X*x2 - 16*X*x3*y1^2 + 32*X*x3*y1*y3 + 16*X*x3*y2^2 - 32*X*x3*y2*y4 - 16*X*x3*y3^2 + 32*X*x3*y3 + 16*X*x3*y4^2 - 8*X*x3 + 16*X*x4*y1^2 - 32*X*x4*y1*y3 - 16*X*x4*y1 - 16*X*x4*y2^2 + 32*X*x4*y2*y4 + 16*X*x4*y2 + 16*X*x4*y3^2 - 16*X*x4*y3 - 16*X*x4*y4^2 - 16*X*x4*y4 - 8*X*x4 + 16*Y^2*x1^2 - 32*Y^2*x1*x2 + 32*Y^2*x1*x3 - 32*Y^2*x1*x4 + 16*Y^2*x2^2 - 32*Y^2*x2*x3 + 32*Y^2*x2*x4 + 16*Y^2*x3^2 - 32*Y^2*x3*x4 + 16*Y^2*x4^2 - 32*Y*x1^2*y3 - 8*Y*x1^2 + 16*Y*x1*x2*y1 - 16*Y*x1*x2*y2 + 48*Y*x1*x2*y3 + 16*Y*x1*x2*y4 - 32*Y*x1*x3*y1 - 32*Y*x1*x3*y3 + 48*Y*x1*x3 + 16*Y*x1*x4*y1 + 16*Y*x1*x4*y2 + 48*Y*x1*x4*y3 - 16*Y*x1*x4*y4 - 16*Y*x2^2*y1 + 16*Y*x2^2*y2 - 16*Y*x2^2*y3 - 16*Y*x2^2*y4 + 8*Y*x2^2 + 48*Y*x2*x3*y1 - 16*Y*x2*x3*y2 + 16*Y*x2*x3*y3 + 16*Y*x2*x3*y4 - 32*Y*x2*x4*y1 - 32*Y*x2*x4*y3 - 48*Y*x2*x4 - 32*Y*x3^2*y1 - 8*Y*x3^2 + 48*Y*x3*x4*y1 + 16*Y*x3*x4*y2 + 16*Y*x3*x4*y3 - 16*Y*x3*x4*y4 - 16*Y*x4^2*y1 - 16*Y*x4^2*y2 - 16*Y*x4^2*y3 + 16*Y*x4^2*y4 + 8*Y*x4^2 + 16*x1^2*y3^2 + 8*x1^2*y3 + x1^2 - 16*x1*x2*y1*y3 - 4*x1*x2*y1 + 16*x1*x2*y2*y3 + 4*x1*x2*y2 - 16*x1*x2*y3^2 - 16*x1*x2*y3*y4 + 4*x1*x2*y3 - 4*x1*x2*y4 + 2*x1*x2 + 16*x1*x3*y1^2 - 24*x1*x3*y1 - 16*x1*x3*y2^2 + 32*x1*x3*y2*y4 + 16*x1*x3*y3^2 - 24*x1*x3*y3 - 16*x1*x3*y4^2 + 2*x1*x3 - 16*x1*x4*y1*y3 - 4*x1*x4*y1 - 16*x1*x4*y2*y3 - 4*x1*x4*y2 - 16*x1*x4*y3^2 + 16*x1*x4*y3*y4 + 4*x1*x4*y3 + 4*x1*x4*y4 + 2*x1*x4 + 4*x2^2*y1^2 - 8*x2^2*y1*y2 + 8*x2^2*y1*y3 + 8*x2^2*y1*y4 - 4*x2^2*y1 + 4*x2^2*y2^2 - 8*x2^2*y2*y3 - 8*x2^2*y2*y4 + 4*x2^2*y2 + 4*x2^2*y3^2 + 8*x2^2*y3*y4 - 4*x2^2*y3 + 4*x2^2*y4^2 - 4*x2^2*y4 + x2^2 - 16*x2*x3*y1^2 + 16*x2*x3*y1*y2 - 16*x2*x3*y1*y3 - 16*x2*x3*y1*y4 + 4*x2*x3*y1 + 4*x2*x3*y2 - 4*x2*x3*y3 - 4*x2*x3*y4 + 2*x2*x3 - 8*x2*x4*y1^2 + 48*x2*x4*y1*y3 + 24*x2*x4*y1 + 8*x2*x4*y2^2 - 16*x2*x4*y2*y4 - 8*x2*x4*y3^2 + 24*x2*x4*y3 + 8*x2*x4*y4^2 + 2*x2*x4 + 16*x3^2*y1^2 + 8*x3^2*y1 + x3^2 - 16*x3*x4*y1^2 - 16*x3*x4*y1*y2 - 16*x3*x4*y1*y3 + 16*x3*x4*y1*y4 + 4*x3*x4*y1 - 4*x3*x4*y2 - 4*x3*x4*y3 + 4*x3*x4*y4 + 2*x3*x4 + 4*x4^2*y1^2 + 8*x4^2*y1*y2 + 8*x4^2*y1*y3 - 8*x4^2*y1*y4 - 4*x4^2*y1 + 4*x4^2*y2^2 + 8*x4^2*y2*y3 - 8*x4^2*y2*y4 - 4*x4^2*y2 + 4*x4^2*y3^2 - 8*x4^2*y3*y4 - 4*x4^2*y3 + 4*x4^2*y4^2 + 4*x4^2*y4 + x4^2)^(1/2) + 2*x1*y1 - 2*x1*y2 + 2*x1*y3 + 2*x3*y1 + 2*x1*y4 - 4*x2*y3 + 2*x3*y2 - 4*x4*y1 + 2*x3*y3 - 2*x3*y4)/(2*x2 + x1*(2*y1 - 2) - 2*x2*y1 + y4*(2*x1 - 2*x2 + 2*x3 - 2*x4) + x3*(2*y1 - 2*y2 - 2*y3 + 2) - x4*(2*y1 - 2*y2 - 2*y3 + 2) - y2*(2*x1 - 2*x2) - y3*(2*x1 - 2*x2))
 xi(2)=  (4*X - x1 - x2 - x3 - x4 - 4*Y*x1 + 4*Y*x2 - 4*Y*x3 + 4*Y*x4 + (16*X^2 - 32*X*Y*x1 + 32*X*Y*x2 - 32*X*Y*x3 + 32*X*Y*x4 - 16*X*x1*y1^2 + 32*X*x1*y1*y3 + 32*X*x1*y1 + 16*X*x1*y2^2 - 32*X*x1*y2*y4 - 16*X*x1*y3^2 + 16*X*x1*y4^2 - 8*X*x1 + 16*X*x2*y1^2 - 32*X*x2*y1*y3 - 16*X*x2*y1 - 16*X*x2*y2^2 + 32*X*x2*y2*y4 - 16*X*x2*y2 + 16*X*x2*y3^2 - 16*X*x2*y3 - 16*X*x2*y4^2 + 16*X*x2*y4 - 8*X*x2 - 16*X*x3*y1^2 + 32*X*x3*y1*y3 + 16*X*x3*y2^2 - 32*X*x3*y2*y4 - 16*X*x3*y3^2 + 32*X*x3*y3 + 16*X*x3*y4^2 - 8*X*x3 + 16*X*x4*y1^2 - 32*X*x4*y1*y3 - 16*X*x4*y1 - 16*X*x4*y2^2 + 32*X*x4*y2*y4 + 16*X*x4*y2 + 16*X*x4*y3^2 - 16*X*x4*y3 - 16*X*x4*y4^2 - 16*X*x4*y4 - 8*X*x4 + 16*Y^2*x1^2 - 32*Y^2*x1*x2 + 32*Y^2*x1*x3 - 32*Y^2*x1*x4 + 16*Y^2*x2^2 - 32*Y^2*x2*x3 + 32*Y^2*x2*x4 + 16*Y^2*x3^2 - 32*Y^2*x3*x4 + 16*Y^2*x4^2 - 32*Y*x1^2*y3 - 8*Y*x1^2 + 16*Y*x1*x2*y1 - 16*Y*x1*x2*y2 + 48*Y*x1*x2*y3 + 16*Y*x1*x2*y4 - 32*Y*x1*x3*y1 - 32*Y*x1*x3*y3 + 48*Y*x1*x3 + 16*Y*x1*x4*y1 + 16*Y*x1*x4*y2 + 48*Y*x1*x4*y3 - 16*Y*x1*x4*y4 - 16*Y*x2^2*y1 + 16*Y*x2^2*y2 - 16*Y*x2^2*y3 - 16*Y*x2^2*y4 + 8*Y*x2^2 + 48*Y*x2*x3*y1 - 16*Y*x2*x3*y2 + 16*Y*x2*x3*y3 + 16*Y*x2*x3*y4 - 32*Y*x2*x4*y1 - 32*Y*x2*x4*y3 - 48*Y*x2*x4 - 32*Y*x3^2*y1 - 8*Y*x3^2 + 48*Y*x3*x4*y1 + 16*Y*x3*x4*y2 + 16*Y*x3*x4*y3 - 16*Y*x3*x4*y4 - 16*Y*x4^2*y1 - 16*Y*x4^2*y2 - 16*Y*x4^2*y3 + 16*Y*x4^2*y4 + 8*Y*x4^2 + 16*x1^2*y3^2 + 8*x1^2*y3 + x1^2 - 16*x1*x2*y1*y3 - 4*x1*x2*y1 + 16*x1*x2*y2*y3 + 4*x1*x2*y2 - 16*x1*x2*y3^2 - 16*x1*x2*y3*y4 + 4*x1*x2*y3 - 4*x1*x2*y4 + 2*x1*x2 + 16*x1*x3*y1^2 - 24*x1*x3*y1 - 16*x1*x3*y2^2 + 32*x1*x3*y2*y4 + 16*x1*x3*y3^2 - 24*x1*x3*y3 - 16*x1*x3*y4^2 + 2*x1*x3 - 16*x1*x4*y1*y3 - 4*x1*x4*y1 - 16*x1*x4*y2*y3 - 4*x1*x4*y2 - 16*x1*x4*y3^2 + 16*x1*x4*y3*y4 + 4*x1*x4*y3 + 4*x1*x4*y4 + 2*x1*x4 + 4*x2^2*y1^2 - 8*x2^2*y1*y2 + 8*x2^2*y1*y3 + 8*x2^2*y1*y4 - 4*x2^2*y1 + 4*x2^2*y2^2 - 8*x2^2*y2*y3 - 8*x2^2*y2*y4 + 4*x2^2*y2 + 4*x2^2*y3^2 + 8*x2^2*y3*y4 - 4*x2^2*y3 + 4*x2^2*y4^2 - 4*x2^2*y4 + x2^2 - 16*x2*x3*y1^2 + 16*x2*x3*y1*y2 - 16*x2*x3*y1*y3 - 16*x2*x3*y1*y4 + 4*x2*x3*y1 + 4*x2*x3*y2 - 4*x2*x3*y3 - 4*x2*x3*y4 + 2*x2*x3 - 8*x2*x4*y1^2 + 48*x2*x4*y1*y3 + 24*x2*x4*y1 + 8*x2*x4*y2^2 - 16*x2*x4*y2*y4 - 8*x2*x4*y3^2 + 24*x2*x4*y3 + 8*x2*x4*y4^2 + 2*x2*x4 + 16*x3^2*y1^2 + 8*x3^2*y1 + x3^2 - 16*x3*x4*y1^2 - 16*x3*x4*y1*y2 - 16*x3*x4*y1*y3 + 16*x3*x4*y1*y4 + 4*x3*x4*y1 - 4*x3*x4*y2 - 4*x3*x4*y3 + 4*x3*x4*y4 + 2*x3*x4 + 4*x4^2*y1^2 + 8*x4^2*y1*y2 + 8*x4^2*y1*y3 - 8*x4^2*y1*y4 - 4*x4^2*y1 + 4*x4^2*y2^2 + 8*x4^2*y2*y3 - 8*x4^2*y2*y4 - 4*x4^2*y2 + 4*x4^2*y3^2 - 8*x4^2*y3*y4 - 4*x4^2*y3 + 4*x4^2*y4^2 + 4*x4^2*y4 + x4^2)^(1/2) + 2*x1*y1 - 2*x1*y2 + 2*x1*y3 + 2*x3*y1 + 2*x1*y4 - 4*x2*y3 + 2*x3*y2 - 4*x4*y1 + 2*x3*y3 - 2*x3*y4)/(2*x2 + x1*(2*y1 - 2) - 2*x2*y1 + y4*(2*x1 - 2*x2 + 2*x3 - 2*x4) + x3*(2*y1 - 2*y2 - 2*y3 + 2) - x4*(2*y1 - 2*y2 - 2*y3 + 2) - y2*(2*x1 - 2*x2) - y3*(2*x1 - 2*x2))

 xi=[...
  -(x1 - 4*X + x2 + x3 + x4 + 4*Y*x1 - 4*Y*x2 + 4*Y*x3 - 4*Y*x4 + (16*X^2 - 32*X*Y*x1 + 32*X*Y*x2 - 32*X*Y*x3 + 32*X*Y*x4 - 16*X*x1*y1^2 + 32*X*x1*y1*y3 + 32*X*x1*y1 + 16*X*x1*y2^2 - 32*X*x1*y2*y4 - 16*X*x1*y3^2 + 16*X*x1*y4^2 - 8*X*x1 + 16*X*x2*y1^2 - 32*X*x2*y1*y3 - 16*X*x2*y1 - 16*X*x2*y2^2 + 32*X*x2*y2*y4 - 16*X*x2*y2 + 16*X*x2*y3^2 - 16*X*x2*y3 - 16*X*x2*y4^2 + 16*X*x2*y4 - 8*X*x2 - 16*X*x3*y1^2 + 32*X*x3*y1*y3 + 16*X*x3*y2^2 - 32*X*x3*y2*y4 - 16*X*x3*y3^2 + 32*X*x3*y3 + 16*X*x3*y4^2 - 8*X*x3 + 16*X*x4*y1^2 - 32*X*x4*y1*y3 - 16*X*x4*y1 - 16*X*x4*y2^2 + 32*X*x4*y2*y4 + 16*X*x4*y2 + 16*X*x4*y3^2 - 16*X*x4*y3 - 16*X*x4*y4^2 - 16*X*x4*y4 - 8*X*x4 + 16*Y^2*x1^2 - 32*Y^2*x1*x2 + 32*Y^2*x1*x3 - 32*Y^2*x1*x4 + 16*Y^2*x2^2 - 32*Y^2*x2*x3 + 32*Y^2*x2*x4 + 16*Y^2*x3^2 - 32*Y^2*x3*x4 + 16*Y^2*x4^2 - 32*Y*x1^2*y3 - 8*Y*x1^2 + 16*Y*x1*x2*y1 - 16*Y*x1*x2*y2 + 48*Y*x1*x2*y3 + 16*Y*x1*x2*y4 - 32*Y*x1*x3*y1 - 32*Y*x1*x3*y3 + 48*Y*x1*x3 + 16*Y*x1*x4*y1 + 16*Y*x1*x4*y2 + 48*Y*x1*x4*y3 - 16*Y*x1*x4*y4 - 16*Y*x2^2*y1 + 16*Y*x2^2*y2 - 16*Y*x2^2*y3 - 16*Y*x2^2*y4 + 8*Y*x2^2 + 48*Y*x2*x3*y1 - 16*Y*x2*x3*y2 + 16*Y*x2*x3*y3 + 16*Y*x2*x3*y4 - 32*Y*x2*x4*y1 - 32*Y*x2*x4*y3 - 48*Y*x2*x4 - 32*Y*x3^2*y1 - 8*Y*x3^2 + 48*Y*x3*x4*y1 + 16*Y*x3*x4*y2 + 16*Y*x3*x4*y3 - 16*Y*x3*x4*y4 - 16*Y*x4^2*y1 - 16*Y*x4^2*y2 - 16*Y*x4^2*y3 + 16*Y*x4^2*y4 + 8*Y*x4^2 + 16*x1^2*y3^2 + 8*x1^2*y3 + x1^2 - 16*x1*x2*y1*y3 - 4*x1*x2*y1 + 16*x1*x2*y2*y3 + 4*x1*x2*y2 - 16*x1*x2*y3^2 - 16*x1*x2*y3*y4 + 4*x1*x2*y3 - 4*x1*x2*y4 + 2*x1*x2 + 16*x1*x3*y1^2 - 24*x1*x3*y1 - 16*x1*x3*y2^2 + 32*x1*x3*y2*y4 + 16*x1*x3*y3^2 - 24*x1*x3*y3 - 16*x1*x3*y4^2 + 2*x1*x3 - 16*x1*x4*y1*y3 - 4*x1*x4*y1 - 16*x1*x4*y2*y3 - 4*x1*x4*y2 - 16*x1*x4*y3^2 + 16*x1*x4*y3*y4 + 4*x1*x4*y3 + 4*x1*x4*y4 + 2*x1*x4 + 4*x2^2*y1^2 - 8*x2^2*y1*y2 + 8*x2^2*y1*y3 + 8*x2^2*y1*y4 - 4*x2^2*y1 + 4*x2^2*y2^2 - 8*x2^2*y2*y3 - 8*x2^2*y2*y4 + 4*x2^2*y2 + 4*x2^2*y3^2 + 8*x2^2*y3*y4 - 4*x2^2*y3 + 4*x2^2*y4^2 - 4*x2^2*y4 + x2^2 - 16*x2*x3*y1^2 + 16*x2*x3*y1*y2 - 16*x2*x3*y1*y3 - 16*x2*x3*y1*y4 + 4*x2*x3*y1 + 4*x2*x3*y2 - 4*x2*x3*y3 - 4*x2*x3*y4 + 2*x2*x3 - 8*x2*x4*y1^2 + 48*x2*x4*y1*y3 + 24*x2*x4*y1 + 8*x2*x4*y2^2 - 16*x2*x4*y2*y4 - 8*x2*x4*y3^2 + 24*x2*x4*y3 + 8*x2*x4*y4^2 + 2*x2*x4 + 16*x3^2*y1^2 + 8*x3^2*y1 + x3^2 - 16*x3*x4*y1^2 - 16*x3*x4*y1*y2 - 16*x3*x4*y1*y3 + 16*x3*x4*y1*y4 + 4*x3*x4*y1 - 4*x3*x4*y2 - 4*x3*x4*y3 + 4*x3*x4*y4 + 2*x3*x4 + 4*x4^2*y1^2 + 8*x4^2*y1*y2 + 8*x4^2*y1*y3 - 8*x4^2*y1*y4 - 4*x4^2*y1 + 4*x4^2*y2^2 + 8*x4^2*y2*y3 - 8*x4^2*y2*y4 - 4*x4^2*y2 + 4*x4^2*y3^2 - 8*x4^2*y3*y4 - 4*x4^2*y3 + 4*x4^2*y4^2 + 4*x4^2*y4 + x4^2)^(1/2) - 2*x1*y1 + 2*x1*y2 - 2*x1*y3 - 2*x3*y1 - 2*x1*y4 + 4*x2*y3 - 2*x3*y2 + 4*x4*y1 - 2*x3*y3 + 2*x3*y4)/(2*x2 + x1*(2*y1 - 2) - 2*x2*y1 + y4*(2*x1 - 2*x2 + 2*x3 - 2*x4) + x3*(2*y1 - 2*y2 - 2*y3 + 2) - x4*(2*y1 - 2*y2 - 2*y3 + 2) - y2*(2*x1 - 2*x2) - y3*(2*x1 - 2*x2))
  (4*X - x1 - x2 - x3 - x4 - 4*Y*x1 + 4*Y*x2 - 4*Y*x3 + 4*Y*x4 + (16*X^2 - 32*X*Y*x1 + 32*X*Y*x2 - 32*X*Y*x3 + 32*X*Y*x4 - 16*X*x1*y1^2 + 32*X*x1*y1*y3 + 32*X*x1*y1 + 16*X*x1*y2^2 - 32*X*x1*y2*y4 - 16*X*x1*y3^2 + 16*X*x1*y4^2 - 8*X*x1 + 16*X*x2*y1^2 - 32*X*x2*y1*y3 - 16*X*x2*y1 - 16*X*x2*y2^2 + 32*X*x2*y2*y4 - 16*X*x2*y2 + 16*X*x2*y3^2 - 16*X*x2*y3 - 16*X*x2*y4^2 + 16*X*x2*y4 - 8*X*x2 - 16*X*x3*y1^2 + 32*X*x3*y1*y3 + 16*X*x3*y2^2 - 32*X*x3*y2*y4 - 16*X*x3*y3^2 + 32*X*x3*y3 + 16*X*x3*y4^2 - 8*X*x3 + 16*X*x4*y1^2 - 32*X*x4*y1*y3 - 16*X*x4*y1 - 16*X*x4*y2^2 + 32*X*x4*y2*y4 + 16*X*x4*y2 + 16*X*x4*y3^2 - 16*X*x4*y3 - 16*X*x4*y4^2 - 16*X*x4*y4 - 8*X*x4 + 16*Y^2*x1^2 - 32*Y^2*x1*x2 + 32*Y^2*x1*x3 - 32*Y^2*x1*x4 + 16*Y^2*x2^2 - 32*Y^2*x2*x3 + 32*Y^2*x2*x4 + 16*Y^2*x3^2 - 32*Y^2*x3*x4 + 16*Y^2*x4^2 - 32*Y*x1^2*y3 - 8*Y*x1^2 + 16*Y*x1*x2*y1 - 16*Y*x1*x2*y2 + 48*Y*x1*x2*y3 + 16*Y*x1*x2*y4 - 32*Y*x1*x3*y1 - 32*Y*x1*x3*y3 + 48*Y*x1*x3 + 16*Y*x1*x4*y1 + 16*Y*x1*x4*y2 + 48*Y*x1*x4*y3 - 16*Y*x1*x4*y4 - 16*Y*x2^2*y1 + 16*Y*x2^2*y2 - 16*Y*x2^2*y3 - 16*Y*x2^2*y4 + 8*Y*x2^2 + 48*Y*x2*x3*y1 - 16*Y*x2*x3*y2 + 16*Y*x2*x3*y3 + 16*Y*x2*x3*y4 - 32*Y*x2*x4*y1 - 32*Y*x2*x4*y3 - 48*Y*x2*x4 - 32*Y*x3^2*y1 - 8*Y*x3^2 + 48*Y*x3*x4*y1 + 16*Y*x3*x4*y2 + 16*Y*x3*x4*y3 - 16*Y*x3*x4*y4 - 16*Y*x4^2*y1 - 16*Y*x4^2*y2 - 16*Y*x4^2*y3 + 16*Y*x4^2*y4 + 8*Y*x4^2 + 16*x1^2*y3^2 + 8*x1^2*y3 + x1^2 - 16*x1*x2*y1*y3 - 4*x1*x2*y1 + 16*x1*x2*y2*y3 + 4*x1*x2*y2 - 16*x1*x2*y3^2 - 16*x1*x2*y3*y4 + 4*x1*x2*y3 - 4*x1*x2*y4 + 2*x1*x2 + 16*x1*x3*y1^2 - 24*x1*x3*y1 - 16*x1*x3*y2^2 + 32*x1*x3*y2*y4 + 16*x1*x3*y3^2 - 24*x1*x3*y3 - 16*x1*x3*y4^2 + 2*x1*x3 - 16*x1*x4*y1*y3 - 4*x1*x4*y1 - 16*x1*x4*y2*y3 - 4*x1*x4*y2 - 16*x1*x4*y3^2 + 16*x1*x4*y3*y4 + 4*x1*x4*y3 + 4*x1*x4*y4 + 2*x1*x4 + 4*x2^2*y1^2 - 8*x2^2*y1*y2 + 8*x2^2*y1*y3 + 8*x2^2*y1*y4 - 4*x2^2*y1 + 4*x2^2*y2^2 - 8*x2^2*y2*y3 - 8*x2^2*y2*y4 + 4*x2^2*y2 + 4*x2^2*y3^2 + 8*x2^2*y3*y4 - 4*x2^2*y3 + 4*x2^2*y4^2 - 4*x2^2*y4 + x2^2 - 16*x2*x3*y1^2 + 16*x2*x3*y1*y2 - 16*x2*x3*y1*y3 - 16*x2*x3*y1*y4 + 4*x2*x3*y1 + 4*x2*x3*y2 - 4*x2*x3*y3 - 4*x2*x3*y4 + 2*x2*x3 - 8*x2*x4*y1^2 + 48*x2*x4*y1*y3 + 24*x2*x4*y1 + 8*x2*x4*y2^2 - 16*x2*x4*y2*y4 - 8*x2*x4*y3^2 + 24*x2*x4*y3 + 8*x2*x4*y4^2 + 2*x2*x4 + 16*x3^2*y1^2 + 8*x3^2*y1 + x3^2 - 16*x3*x4*y1^2 - 16*x3*x4*y1*y2 - 16*x3*x4*y1*y3 + 16*x3*x4*y1*y4 + 4*x3*x4*y1 - 4*x3*x4*y2 - 4*x3*x4*y3 + 4*x3*x4*y4 + 2*x3*x4 + 4*x4^2*y1^2 + 8*x4^2*y1*y2 + 8*x4^2*y1*y3 - 8*x4^2*y1*y4 - 4*x4^2*y1 + 4*x4^2*y2^2 + 8*x4^2*y2*y3 - 8*x4^2*y2*y4 - 4*x4^2*y2 + 4*x4^2*y3^2 - 8*x4^2*y3*y4 - 4*x4^2*y3 + 4*x4^2*y4^2 + 4*x4^2*y4 + x4^2)^(1/2) + 2*x1*y1 - 2*x1*y2 + 2*x1*y3 + 2*x3*y1 + 2*x1*y4 - 4*x2*y3 + 2*x3*y2 - 4*x4*y1 + 2*x3*y3 - 2*x3*y4)/(2*x2 + x1*(2*y1 - 2) - 2*x2*y1 + y4*(2*x1 - 2*x2 + 2*x3 - 2*x4) + x3*(2*y1 - 2*y2 - 2*y3 + 2) - x4*(2*y1 - 2*y2 - 2*y3 + 2) - y2*(2*x1 - 2*x2) - y3*(2*x1 - 2*x2))
 -(4*X - x1 - x2 - x3 - x4 - 4*Y*x1 + 4*Y*x2 - 4*Y*x3 + 4*Y*x4 + (16*X^2 - 32*X*Y*x1 + 32*X*Y*x2 - 32*X*Y*x3 + 32*X*Y*x4 - 16*X*x1*y1^2 + 32*X*x1*y1*y3 + 32*X*x1*y1 + 16*X*x1*y2^2 - 32*X*x1*y2*y4 - 16*X*x1*y3^2 + 16*X*x1*y4^2 - 8*X*x1 + 16*X*x2*y1^2 - 32*X*x2*y1*y3 - 16*X*x2*y1 - 16*X*x2*y2^2 + 32*X*x2*y2*y4 - 16*X*x2*y2 + 16*X*x2*y3^2 - 16*X*x2*y3 - 16*X*x2*y4^2 + 16*X*x2*y4 - 8*X*x2 - 16*X*x3*y1^2 + 32*X*x3*y1*y3 + 16*X*x3*y2^2 - 32*X*x3*y2*y4 - 16*X*x3*y3^2 + 32*X*x3*y3 + 16*X*x3*y4^2 - 8*X*x3 + 16*X*x4*y1^2 - 32*X*x4*y1*y3 - 16*X*x4*y1 - 16*X*x4*y2^2 + 32*X*x4*y2*y4 + 16*X*x4*y2 + 16*X*x4*y3^2 - 16*X*x4*y3 - 16*X*x4*y4^2 - 16*X*x4*y4 - 8*X*x4 + 16*Y^2*x1^2 - 32*Y^2*x1*x2 + 32*Y^2*x1*x3 - 32*Y^2*x1*x4 + 16*Y^2*x2^2 - 32*Y^2*x2*x3 + 32*Y^2*x2*x4 + 16*Y^2*x3^2 - 32*Y^2*x3*x4 + 16*Y^2*x4^2 - 32*Y*x1^2*y3 - 8*Y*x1^2 + 16*Y*x1*x2*y1 - 16*Y*x1*x2*y2 + 48*Y*x1*x2*y3 + 16*Y*x1*x2*y4 - 32*Y*x1*x3*y1 - 32*Y*x1*x3*y3 + 48*Y*x1*x3 + 16*Y*x1*x4*y1 + 16*Y*x1*x4*y2 + 48*Y*x1*x4*y3 - 16*Y*x1*x4*y4 - 16*Y*x2^2*y1 + 16*Y*x2^2*y2 - 16*Y*x2^2*y3 - 16*Y*x2^2*y4 + 8*Y*x2^2 + 48*Y*x2*x3*y1 - 16*Y*x2*x3*y2 + 16*Y*x2*x3*y3 + 16*Y*x2*x3*y4 - 32*Y*x2*x4*y1 - 32*Y*x2*x4*y3 - 48*Y*x2*x4 - 32*Y*x3^2*y1 - 8*Y*x3^2 + 48*Y*x3*x4*y1 + 16*Y*x3*x4*y2 + 16*Y*x3*x4*y3 - 16*Y*x3*x4*y4 - 16*Y*x4^2*y1 - 16*Y*x4^2*y2 - 16*Y*x4^2*y3 + 16*Y*x4^2*y4 + 8*Y*x4^2 + 16*x1^2*y3^2 + 8*x1^2*y3 + x1^2 - 16*x1*x2*y1*y3 - 4*x1*x2*y1 + 16*x1*x2*y2*y3 + 4*x1*x2*y2 - 16*x1*x2*y3^2 - 16*x1*x2*y3*y4 + 4*x1*x2*y3 - 4*x1*x2*y4 + 2*x1*x2 + 16*x1*x3*y1^2 - 24*x1*x3*y1 - 16*x1*x3*y2^2 + 32*x1*x3*y2*y4 + 16*x1*x3*y3^2 - 24*x1*x3*y3 - 16*x1*x3*y4^2 + 2*x1*x3 - 16*x1*x4*y1*y3 - 4*x1*x4*y1 - 16*x1*x4*y2*y3 - 4*x1*x4*y2 - 16*x1*x4*y3^2 + 16*x1*x4*y3*y4 + 4*x1*x4*y3 + 4*x1*x4*y4 + 2*x1*x4 + 4*x2^2*y1^2 - 8*x2^2*y1*y2 + 8*x2^2*y1*y3 + 8*x2^2*y1*y4 - 4*x2^2*y1 + 4*x2^2*y2^2 - 8*x2^2*y2*y3 - 8*x2^2*y2*y4 + 4*x2^2*y2 + 4*x2^2*y3^2 + 8*x2^2*y3*y4 - 4*x2^2*y3 + 4*x2^2*y4^2 - 4*x2^2*y4 + x2^2 - 16*x2*x3*y1^2 + 16*x2*x3*y1*y2 - 16*x2*x3*y1*y3 - 16*x2*x3*y1*y4 + 4*x2*x3*y1 + 4*x2*x3*y2 - 4*x2*x3*y3 - 4*x2*x3*y4 + 2*x2*x3 - 8*x2*x4*y1^2 + 48*x2*x4*y1*y3 + 24*x2*x4*y1 + 8*x2*x4*y2^2 - 16*x2*x4*y2*y4 - 8*x2*x4*y3^2 + 24*x2*x4*y3 + 8*x2*x4*y4^2 + 2*x2*x4 + 16*x3^2*y1^2 + 8*x3^2*y1 + x3^2 - 16*x3*x4*y1^2 - 16*x3*x4*y1*y2 - 16*x3*x4*y1*y3 + 16*x3*x4*y1*y4 + 4*x3*x4*y1 - 4*x3*x4*y2 - 4*x3*x4*y3 + 4*x3*x4*y4 + 2*x3*x4 + 4*x4^2*y1^2 + 8*x4^2*y1*y2 + 8*x4^2*y1*y3 - 8*x4^2*y1*y4 - 4*x4^2*y1 + 4*x4^2*y2^2 + 8*x4^2*y2*y3 - 8*x4^2*y2*y4 - 4*x4^2*y2 + 4*x4^2*y3^2 - 8*x4^2*y3*y4 - 4*x4^2*y3 + 4*x4^2*y4^2 + 4*x4^2*y4 + x4^2)^(1/2) + 2*x1*y1 + 2*x1*y2 - 4*x2*y1 + 2*x1*y3 + 2*x3*y1 - 2*x1*y4 - 2*x3*y2 + 2*x3*y3 + 2*x3*y4 - 4*x4*y3)/(2*(x1 + x2 - x3 - x4 - x1*y1 - x1*y2 + x2*y1 + x1*y3 + x2*y2 - x3*y1 + x1*y4 - x2*y3 - x3*y2 + x4*y1 - x2*y4 + x3*y3 + x4*y2 + x3*y4 - x4*y3 - x4*y4))
  (x1 - 4*X + x2 + x3 + x4 + 4*Y*x1 - 4*Y*x2 + 4*Y*x3 - 4*Y*x4 + (16*X^2 - 32*X*Y*x1 + 32*X*Y*x2 - 32*X*Y*x3 + 32*X*Y*x4 - 16*X*x1*y1^2 + 32*X*x1*y1*y3 + 32*X*x1*y1 + 16*X*x1*y2^2 - 32*X*x1*y2*y4 - 16*X*x1*y3^2 + 16*X*x1*y4^2 - 8*X*x1 + 16*X*x2*y1^2 - 32*X*x2*y1*y3 - 16*X*x2*y1 - 16*X*x2*y2^2 + 32*X*x2*y2*y4 - 16*X*x2*y2 + 16*X*x2*y3^2 - 16*X*x2*y3 - 16*X*x2*y4^2 + 16*X*x2*y4 - 8*X*x2 - 16*X*x3*y1^2 + 32*X*x3*y1*y3 + 16*X*x3*y2^2 - 32*X*x3*y2*y4 - 16*X*x3*y3^2 + 32*X*x3*y3 + 16*X*x3*y4^2 - 8*X*x3 + 16*X*x4*y1^2 - 32*X*x4*y1*y3 - 16*X*x4*y1 - 16*X*x4*y2^2 + 32*X*x4*y2*y4 + 16*X*x4*y2 + 16*X*x4*y3^2 - 16*X*x4*y3 - 16*X*x4*y4^2 - 16*X*x4*y4 - 8*X*x4 + 16*Y^2*x1^2 - 32*Y^2*x1*x2 + 32*Y^2*x1*x3 - 32*Y^2*x1*x4 + 16*Y^2*x2^2 - 32*Y^2*x2*x3 + 32*Y^2*x2*x4 + 16*Y^2*x3^2 - 32*Y^2*x3*x4 + 16*Y^2*x4^2 - 32*Y*x1^2*y3 - 8*Y*x1^2 + 16*Y*x1*x2*y1 - 16*Y*x1*x2*y2 + 48*Y*x1*x2*y3 + 16*Y*x1*x2*y4 - 32*Y*x1*x3*y1 - 32*Y*x1*x3*y3 + 48*Y*x1*x3 + 16*Y*x1*x4*y1 + 16*Y*x1*x4*y2 + 48*Y*x1*x4*y3 - 16*Y*x1*x4*y4 - 16*Y*x2^2*y1 + 16*Y*x2^2*y2 - 16*Y*x2^2*y3 - 16*Y*x2^2*y4 + 8*Y*x2^2 + 48*Y*x2*x3*y1 - 16*Y*x2*x3*y2 + 16*Y*x2*x3*y3 + 16*Y*x2*x3*y4 - 32*Y*x2*x4*y1 - 32*Y*x2*x4*y3 - 48*Y*x2*x4 - 32*Y*x3^2*y1 - 8*Y*x3^2 + 48*Y*x3*x4*y1 + 16*Y*x3*x4*y2 + 16*Y*x3*x4*y3 - 16*Y*x3*x4*y4 - 16*Y*x4^2*y1 - 16*Y*x4^2*y2 - 16*Y*x4^2*y3 + 16*Y*x4^2*y4 + 8*Y*x4^2 + 16*x1^2*y3^2 + 8*x1^2*y3 + x1^2 - 16*x1*x2*y1*y3 - 4*x1*x2*y1 + 16*x1*x2*y2*y3 + 4*x1*x2*y2 - 16*x1*x2*y3^2 - 16*x1*x2*y3*y4 + 4*x1*x2*y3 - 4*x1*x2*y4 + 2*x1*x2 + 16*x1*x3*y1^2 - 24*x1*x3*y1 - 16*x1*x3*y2^2 + 32*x1*x3*y2*y4 + 16*x1*x3*y3^2 - 24*x1*x3*y3 - 16*x1*x3*y4^2 + 2*x1*x3 - 16*x1*x4*y1*y3 - 4*x1*x4*y1 - 16*x1*x4*y2*y3 - 4*x1*x4*y2 - 16*x1*x4*y3^2 + 16*x1*x4*y3*y4 + 4*x1*x4*y3 + 4*x1*x4*y4 + 2*x1*x4 + 4*x2^2*y1^2 - 8*x2^2*y1*y2 + 8*x2^2*y1*y3 + 8*x2^2*y1*y4 - 4*x2^2*y1 + 4*x2^2*y2^2 - 8*x2^2*y2*y3 - 8*x2^2*y2*y4 + 4*x2^2*y2 + 4*x2^2*y3^2 + 8*x2^2*y3*y4 - 4*x2^2*y3 + 4*x2^2*y4^2 - 4*x2^2*y4 + x2^2 - 16*x2*x3*y1^2 + 16*x2*x3*y1*y2 - 16*x2*x3*y1*y3 - 16*x2*x3*y1*y4 + 4*x2*x3*y1 + 4*x2*x3*y2 - 4*x2*x3*y3 - 4*x2*x3*y4 + 2*x2*x3 - 8*x2*x4*y1^2 + 48*x2*x4*y1*y3 + 24*x2*x4*y1 + 8*x2*x4*y2^2 - 16*x2*x4*y2*y4 - 8*x2*x4*y3^2 + 24*x2*x4*y3 + 8*x2*x4*y4^2 + 2*x2*x4 + 16*x3^2*y1^2 + 8*x3^2*y1 + x3^2 - 16*x3*x4*y1^2 - 16*x3*x4*y1*y2 - 16*x3*x4*y1*y3 + 16*x3*x4*y1*y4 + 4*x3*x4*y1 - 4*x3*x4*y2 - 4*x3*x4*y3 + 4*x3*x4*y4 + 2*x3*x4 + 4*x4^2*y1^2 + 8*x4^2*y1*y2 + 8*x4^2*y1*y3 - 8*x4^2*y1*y4 - 4*x4^2*y1 + 4*x4^2*y2^2 + 8*x4^2*y2*y3 - 8*x4^2*y2*y4 - 4*x4^2*y2 + 4*x4^2*y3^2 - 8*x4^2*y3*y4 - 4*x4^2*y3 + 4*x4^2*y4^2 + 4*x4^2*y4 + x4^2)^(1/2) - 2*x1*y1 - 2*x1*y2 + 4*x2*y1 - 2*x1*y3 - 2*x3*y1 + 2*x1*y4 + 2*x3*y2 - 2*x3*y3 - 2*x3*y4 + 4*x4*y3)/(2*(x1 + x2 - x3 - x4 - x1*y1 - x1*y2 + x2*y1 + x1*y3 + x2*y2 - x3*y1 + x1*y4 - x2*y3 - x3*y2 + x4*y1 - x2*y4 + x3*y3 + x4*y2 + x3*y4 - x4*y3 - x4*y4))];
 
end