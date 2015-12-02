function [D,M]=gencouplingmatrices(params,elements,nodes)

imele=find(params.ifaceele==2);
isele=find(params.ifaceele==1);

imnode=find(params.ifacenode==2);
isnode=find(params.ifacenode==1);

numnodes   =length(nodes);
numdof     =numnodes*2;
numlagnodes=length(isnode);
numlagdof  =numlagnodes*2;

nummodnodes=length(params.modnodes)
nummoddof=2*nummodnodes;

[gpxi,gpw] =gaussrulequad();
numgp=length(gpw)



% 'please check this for correctness'
pause(1)

D=zeros(numdof+numlagdof,...
        numdof+numlagdof);
M=zeros(numdof+numlagdof,...
        numdof+numlagdof);

      
D=zeros(numdof+numlagdof+nummoddof,...
        numdof+numlagdof+nummoddof);
M=zeros(numdof+numlagdof+nummoddof,...
        numdof+numlagdof+nummoddof);

%M matrix
for il=isele
    
    nodeslag=nodes(elements(il,:),:);
    dofslag=dofs(elements(il,:))+params.dofshift;
    curlagele=elements(il,:);
    
    
    
    for im=imele
        curmasterele=elements(im,:);
        dofsmaster=dofs(curmasterele);
        nodemaster=nodes(elements(im,:),:);
        
        Me=zeros(8,8);
        for igp=1:numgp
            xilag=gpxi(igp,:);
            Xglob=LocalToGlobal( xilag,nodeslag );
            ximaster=GlobalToLocal( Xglob,nodemaster );
            if CheckMapping('quad4',ximaster)==1
                Nl =sval('quad4',xilag);
                

                if sum(params.modele==il)>0 %% this element uses modified 
                    Nl=svalmod(elements(il,:),params,xilag);
                end
                
                Nm =sval('quad4',ximaster);

                Nlmat=[...
                      Nl(1) 0     Nl(2) 0     Nl(3) 0     Nl(4) 0
                      0     Nl(1) 0     Nl(2) 0     Nl(3) 0     Nl(4)];
                Nmmat=[...
                      Nm(1) 0     Nm(2) 0     Nm(3) 0     Nm(4) 0
                      0     Nm(1) 0     Nm(2) 0     Nm(3) 0     Nm(4)];
                    
                jac=det(Jacobian(nodeslag,xilag));
                Me=Me+Nlmat'*Nmmat*gpw(igp)*inv(jac); %only valid because I am using matching interfaces
            end
        end
        M(dofslag,dofsmaster)=M(dofslag,dofsmaster)+Me;
    end
    
end



%D matrix
for il=isele
    
    nodelag=nodes(elements(il,:),:);
    dofslag=dofs(elements(il,:))+params.dofshift;
    curlagele=elements(il,:);
    
    
    
    for is=isele
        curslave=elements(is,:);
        dofsslave=dofs(curslave);
        nodeslave=nodes(elements(is,:),:);
        
        De=zeros(8,8);
        for igp=1:numgp
            xilag=gpxi(igp,:);
            Xglob=LocalToGlobal( xilag,nodelag );
            xislave=GlobalToLocal( Xglob,nodeslave );
            if CheckMapping('quad4',xislave)==1
                Nl =sval('quad4',xilag);
                T=[4 -2  1 -2
                  -2  4 -2  1
                   1 -2  4 -2
                  -2  1 -2  4];
                %Nl=T*Nl;
                
                if sum(params.modele==il)>0 %% this element uses modified 
                    Nl=svalmod(elements(il,:),params,xilag);
                end
                
                Ns =sval('quad4',xislave);

                Nlmat=[...
                      Nl(1) 0     Nl(2) 0     Nl(3) 0     Nl(4) 0
                      0     Nl(1) 0     Nl(2) 0     Nl(3) 0     Nl(4)];
                Nsmat=[...
                      Ns(1) 0     Ns(2) 0     Ns(3) 0     Ns(4) 0
                      0     Ns(1) 0     Ns(2) 0     Ns(3) 0     Ns(4)];
                jac=det(Jacobian(nodeslag,xilag));
                De=De+Nlmat'*Nsmat*gpw(igp)*inv(jac); %only valid because I am using matching interfaces
            end
        end
        D(dofslag,dofsslave)=D(dofslag,dofsslave)+De;
    end
    
end






% %D matrix
% for is1=isele
%     dofslag=dofs(elements(is1,:))+params.dofshift
%     for i2=isele
%         dofsslave=dofs(elements(is1,:))
%         
%         De=zeros(8,8);
%         for igp=1:numgp
%             N =sval('quad4',gpxi(igp,:))
%             Nmat=[...
%                   N(1)  0     N(2)  0     N(3)  0     N(4)  0
%                   0     N(1)  0     N(2)  0     N(3)  0     N(4)];
%             De=De+Nmat'*Nmat*gpw(igp) %only valid because I am using matching interfaces
%         end
%         D(dofslag,dofsslave)=D(dofslag,dofsslave)+De;
%     end
% end


D=D(1:numdof+numlagdof,1:numdof+numlagdof);
M=M(1:numdof+numlagdof,1:numdof+numlagdof);



end