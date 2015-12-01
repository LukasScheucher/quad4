clear
close all
clc


syms a b c

solve(a+b-c,a-b+3*c+2,[a,b])



syms x1 x2 x3 x4
syms y1 y2 y3 y4

syms X Y


syms r s




sol=solve(...
    (x1+x2+x3+x4)+r*(-x1+x2+x3-x4)+s*(-x1-x2+x3+x4)+r*s*(x1-x2+x3-x4)-4*X,...
    (y1+y2+y3+y4)+r*(-y1+y2+y3-y4)+s*(-y1-y2+y3+y4)+r*s+(y1-y2+y3-y4)-4*Y,[r,s])


rsol=simplify(sol.r)
ssol=simplify(sol.s)
    