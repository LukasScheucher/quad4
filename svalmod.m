function Nmod =svalmod(curnodes,params,xi)


N=sval('quad4',xi);

curmod=[0 0 0 0];

for testnode=params.modnodes
      curmod=curmod+[curnodes==testnode];
end


curmod;
fac=1/sum(curmod==0);
summod=transpose(N)*curmod';

Nmod=N.*(1-curmod');
Nmod=Nmod+(fac*summod).*(1-curmod');



end