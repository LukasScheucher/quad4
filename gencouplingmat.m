function [Cg]=gencouplingmat(params,nodes,lagele,counterele,modlag,modcounter,valids,validm)

   
    TOL=1e-12;
    
    nodeslag=nodes(lagele,:);
    dofslag=dofs(lagele)+params.dofshift;
    
    dofscounter=dofs(counterele);
    nodescounter=nodes(counterele,:);
    
    
    [gpxi,gpw] =gaussrulequad();
    numgp=length(gpw);
    
    
    syssize_ex=length(nodes)*2+2*sum(params.bodynode==1);
    syssize=   length(nodes)*2+2*sum(params.ifacenode==1);
        
       
    curmod=[0 0 0 0]
    for testnode=params.modnodes
      curmod=curmod+[lagele==testnode]
    end

    C=zeros(8,8);
    for igp=1:numgp
      
        xilag=gpxi(igp,:);
        Xglob=LocalToGlobal( xilag,nodeslag );
        xicounter=GlobalToLocal( Xglob,nodescounter );
        if CheckMapping('quad4',xicounter)==1
            
          Nl =sval('quad4',xilag);

%% temporarily disable the modified shit at all
%           if length(intersect(lagele,params.modnodes))>0
%           %%if sum(params.modele==il)>0 %% this element uses modified 
%               Nl=svalmod(elements(il,:),params,xilag);
%           end
%%
           if strcmp(modlag,'mod')
             Nl=[...
                    0
                    0
                    ( 1+xilag(1) )
                    ( 1-xilag(1) )]*0.5;         
           end
          
%           T=[4 -2  1 -2
%           -2  4 -2  1
%            1 -2  4 -2
%           -2  1 -2  4];
%           Nl=T*Nl;

          Nc =sval('quad4',xicounter);
          
          
          if strcmp(modcounter,'mod')
             Nl=[...
                    0
                    0
                    ( 1+xilag(1) )
                    ( 1-xilag(1) )]*0.5;         
           end

          Nlmat=[...
                Nl(1) 0     Nl(2) 0     Nl(3) 0     Nl(4) 0
                0     Nl(1) 0     Nl(2) 0     Nl(3) 0     Nl(4)];
          Nsmat=[...
                Nc(1) 0     Nc(2) 0     Nc(3) 0     Nc(4) 0
                0     Nc(1) 0     Nc(2) 0     Nc(3) 0     Nc(4)];
          jac=det(Jacobian(nodeslag,xilag));
          C=C+Nlmat'*Nsmat*gpw(igp)*inv(jac);
        end
    end
    
    
    C=C.*(C>TOL);
    
    Cg=zeros(syssize_ex,syssize_ex);%extend
    Cg(dofs(lagele)+params.dofshift,dofs(counterele) )=C;
    Cg=Cg(1:syssize,1:syssize);%reduce
   

end