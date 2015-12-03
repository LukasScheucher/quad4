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



%%HACK

% r=xi(1);
% s=xi(2);
% Nmod=[...
%       0
%       0
%       (1+r)
%       (1-r)]*0.5;




%%%% END HACK


%%HACK

% r=xi(1);
% s=xi(2);
% Nmod=[...
%       (1-r)*(1-s)
%       (1+r)*(1-s)
%       (1+r)*(1+s)
%       (1-r)*(1+s)]*0.25;
%%%% END HACK


end