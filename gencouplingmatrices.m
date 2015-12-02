function [D,M]=gencouplingmatrices(params,elements,nodes)

imele=[find(params.ifaceele==2),find(params.ifaceele==3)]; %indices of all master elements
isele=[find(params.ifaceele==1),find(params.ifaceele==3)]; %indices of all slave elements
ilele=[find(params.ifaceele==1),find(params.ifaceele==3)]; %indices of all lagrange elements

imnode=[find(params.ifacenode==2),find(params.ifacenode==3)]; %indices of all master nodes
isnode=[find(params.ifacenode==1),find(params.ifacenode==3)]; %indices of all slave elements
ilnode= find(params.ifacenode==1);
%isnode=[find(params.ifacenode==1),find(params.ifacenode==3)];

numnodes   =length(nodes);
numdof     =numnodes*2;

numlagnodes=length(ilnode);
numlagdof  =numlagnodes*2;


numslavenodes=length(isnode);
numslavedof  =numslavenodes*2;

nummasternodes=length(imnode);
nummasterdof  =nummasternodes*2;


nummodnodes=length(params.modnodes);
nummoddof=2*nummodnodes;



ibmele=find(params.bodyele==2); %indices of all master elements
ibsele=find(params.bodyele==1); %indices of all slave elements




%pause(100)



[gpxi,gpw] =gaussrulequad();
numgp=length(gpw);



D=zeros(numdof+numlagdof,...
        numdof+numlagdof);
M=zeros(numdof+numlagdof,...
        numdof+numlagdof);
      
      
D=zeros(numdof+numslavedof,...
        numdof+numslavedof);
M=zeros(numdof+numslavedof,...
        numdof+numslavedof);


%M matrix
for il=ilele
    
    nodeslag=nodes(elements(il,:),:);
    dofslag=dofs(elements(il,:))+params.dofshift;
    curlagele=elements(il,:);
    
    
    
    for im=ibmele
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now perform the extra coupling of additional M



%M matrix
for il=ilele
    
    nodeslag=nodes(elements(il,:),:);
    dofslag=dofs(elements(il,:))+params.dofshift;
    curlagele=elements(il,:);
    
    
    
    for im=params.modele
        curmasterele=elements(im,:);
        dofsmaster=dofs(curmasterele);
        nodemaster=nodes(elements(im,:),:);
        
        curmasterele
        curmod=[0 0 0 0]
        for testnode=params.modnodes
          curmod=curmod+[curmasterele==testnode]
        end
        
        
        %corrector is ma matrix of only ones and zeros
        % its sole purpose is to delete all
        corrector=ones(8,8);
        for i=1:4
          if curmod(i)==0
            corrector(i*2,:)=0;
            corrector(:,i*2)=0;
            corrector(i*2-1,:)=0;
            corrector(:,i*2-1)=0;
          end
        end
        
        corrector
        
        
        
        %pause(100);
        
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
                
                Me=Me.*corrector;
            end
        end
        M(dofslag,dofsmaster)=M(dofslag,dofsmaster)+Me;
    end
    
end


%% done additional part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%D matrix
for il=ilele
    
    nodelag=nodes(elements(il,:),:);
    dofslag=dofs(elements(il,:))+params.dofshift;
    curlagele=elements(il,:);
    
    
    
    for is=ibsele
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



% Die zusaetyliochen Nullzeilen wurde nur eingefuehrt, um den einfachen
% Assemblierungsprozess beibehalten zu koennen. es sollten nur Nullzeilen
% assembliert worden sein.

if sum(sum(D(numdof+numlagdof:end,numdof+numlagdof:end)+...
           M(numdof+numlagdof:end,numdof+numlagdof:end)))>0
         error('auxilary lines of D and M are not empty');
end
    
D=D(1:numdof+numlagdof,1:numdof+numlagdof);
M=M(1:numdof+numlagdof,1:numdof+numlagdof);



end