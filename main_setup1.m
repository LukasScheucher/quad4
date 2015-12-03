clear
close all
clc

%% setup
setup1
f=figure()
%disp=zeros(length(nodes),2);
plotframe(f,params,elements,nodes)

%% stiffmate
stiff = AssembleStiff(params,elements,nodes);

%% standard coupling matrices
[M73]=gencouplingmat(params,nodes,elements(7,:),elements(3,:),1,1)
[M84]=gencouplingmat(params,nodes,elements(8,:),elements(4,:),1,1)

[D77]=gencouplingmat(params,nodes,elements(7,:),elements(7,:),1,1)
[D88]=gencouplingmat(params,nodes,elements(8,:),elements(8,:),1,1)


%% modified coupling matrices

[M95]=gencouplingmat(params,nodes,elements(9,:),elements(5,:),1,1)
[M106]=gencouplingmat(params,nodes,elements(10,:),elements(6,:),1,1)

[D99]=gencouplingmat(params,nodes,elements(9,:),elements(9,:),1,1)
[D1010]=gencouplingmat(params,nodes,elements(10,:),elements(10,:),1,1)

D=D77+D88+D99+D1010;
M=M73+M84+M95+M106;

%% solution
LHS=stiff+(M-D+M'-D')*1000;
RHS=zeros(size(LHS,1),1) ;

imagesc(LHS==0)
[LHS,RHS]=ApplyDirich(LHS,RHS,dirichdof);

sol=LHS\RHS;

disp=[sol(1:2:numnodes*2),sol(2:2:numnodes*2)];



f=figure();
axis equal
%disp=zeros(21,2)
plotframe(f,params,elements,nodes)

hold on

plotdisp(f,params,elements,nodes,disp*4)

%%
f=figure()
plotframe(f,params,elements,nodes)

hold on
for iele=7:length(elements)

 plotsf(f,elements,nodes,iele,params)


end