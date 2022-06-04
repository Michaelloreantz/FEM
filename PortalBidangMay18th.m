%% Plane Frame Stiffness Method
% Michael Loreantz Steven Tambunan // 1906301923
%% Prerequisites
clear
close all
clc

disp('<strong>Michael Loreantz Steven Tambunan // 1906301923</strong>')
disp('<strong>Plane Frame Stiffness Method</strong>')
%% Data Inputs
DOF = 3;                %Number of DOFs per node
Node = [1 0 0
        2 0 1
        3 0 2
        4 0 3
        5 0 4
        6 1 4
        7 2 4
        8 3 4
        9 4 4
        10 5 4
        11 6 4
        12 7 4
        13 8 4
        14 8.5 3
        15 9 2
        16 9.5 1
        17 10 0
        18 0 5
        19 0 6
        20 0 7
        21 0 8];         %Node locations (nodenumber,x, y) (m)
Member = [1 1 2 1
          2 2 3 1
          3 3 4 1
          4 4 5 1
          5 5 6 1
          6 6 7 1
          7 7 8 1
          8 8 9 1
          9 9 10 1
          10 10 11 1
          11 11 12 1
          12 12 13 1
          13 13 14 1
          14 14 15 1
          15 15 16 1
          16 16 17 1
          17 5 18 1
          18 18 19 1
          19 19 20 1
          20 20 21 1];    %Connectivity (Member Number, i node, j node, 1=Frame, 2=Truss)
Ex = 200*10^6;            %Given Variable E
Ax = 6000*10^-6;          %Given Variable A
Ix = 3*10^-6;             %Given Variable I
E = Ex;                   %Modulus of elasticity ([MemberNumber, Elasticity] OR 1 value that applies for all member) (kN/m^2) 
I = Ix;                   %Moment of inertia along flexure axis ([MemberNumber, Moment of Inertia] OR 1 value that applies for all member) (m^4)
A = Ax;                   %Cross sectional area ([MemberNumber, Sectional Area] OR 1 value that applies for all member) (m^2)
Px=0;                     %Given Variable P
Mx=0;                     %Given Variable M
f0 = 0;                   %Given Variable f0
nDOF = [1 1 1 1
        17 0 1 0
        21 1 1 0];        %Nodes with 1 or more boundary conditions (node, u, v, theta) (non zero value=restrained). Note: For Truss Members, theta=restrained.
    
interloads = [13 40 40 1
              14 40 40 1
              15 40 40 1
              16 40 40 1]; %Intermediate loadings (member, i node intermediate load, j node intermediate load, direction (1=local -y, 2=global -y, 3=global +x))
    
P = [5 0 0 -75];                %External nodal loads (node, fy, M)
FigScale=0.5;                   %Figure Scaling
ShearScale=1/200;               %Shear Diagram Scaling
BendingMomentScale=1/200;       %Bending Moment Diagram Scaling
AxialScale=1/200;               %Axial Diagram Scaling
%% Force vector assembly (F)
% Equivalent nodal loads assembly
memberid=Member(:,1);                                 %List of available members
if isempty(interloads)==1                             %If there is no intermediate loadings then:
    intloadmember=[];                                 %Members with intermediate loading is none
    intloadmemberx=setxor(intloadmember,Member(:,1)); %Members without intermediate loadings is all member
else
    intloadmember=interloads(:,1);                    %Members with intermediate loadings
	intloadmemberx=setxor(intloadmember,Member(:,1)); %Members without intermediate loadings
end
memberid(intloadmemberx)=0;                       %Excluding members without intermediate loadings

elemENL=zeros(size(Member,1),1+2*DOF);            %Zero matrix to save all member equivalent nodal loadings
for elnum=1:size(elemENL,1)                       %Looping for each member
    if memberid(elnum)==elnum                     %If a member contains intermediate loads, then calculate equivalent nodal load
       id2el=find(intloadmember==elnum);          %To help indexing the members with intermediate loading, as not all member will have one.
       inode=Member(elnum,2);                     %i node
       ix=Node(inode,2);                          %i node x coordinate
       iy=Node(inode,3);                          %i node y coordinate
       jnode=Member(elnum,3);                     %j node
       jx=Node(jnode,2);                          %j node x coordinate
       jy=Node(jnode,3);                          %j node y coordinate
       l=sqrt((jx-ix)^2+(jy-iy)^2);               %Member length
       c=(jx-ix)/l;                               %Cosine of element
       s=(jy-iy)/l;   
       
       if interloads(id2el,4)==1                  %If intermediate load is applied downward to member's local y axis, then:
       w1t=interloads(id2el,2);                   %Transverse w1=w1
       w2t=interloads(id2el,3);                   %Transverse w2=w2
       w1a=0;                                     %No axial w1
       w2a=0;                                     %No axial w2
       elseif interloads(id2el,4)==2              %Else if intermediate load is applied downward to global y axis, then:
       w1t=interloads(id2el,2)*c;                 %Transverse w1=w1*cos(theta)
       w2t=interloads(id2el,3)*c;                 %Transverse w2=w2*cos(theta)
       w1a=interloads(id2el,2)*s;                 %Axial w1=w1*sin(theta)
       w2a=interloads(id2el,3)*s;                 %Axial w2=w2*sin(theta)
       elseif interloads(id2el,4)==3              %Else if intermediate load is applied to positive global x axis, then:
       w1t=interloads(id2el,2)*s;                 %Transverse w1=w1*sin(theta)
       w2t=interloads(id2el,3)*s;                 %Transverse w2=w2*sin(theta)
       w1a=-interloads(id2el,2)*c;                %Axial w1=-w1*cos(theta)
       w2a=-interloads(id2el,3)*c;                %Axial w2=-w2*cos(theta)
       end
       [enli,T,enl]=interload(w1t,w2t,w1a,w2a,l,c,s);    %Equivalent nodal loading calculation
       ENL=[elnum enl'];                                %Indexed form of nodal loading calculation
    else
       ENL=[elnum zeros(1,2*DOF)];                      %If a member does not contain intermediate loading, then do not calculate equivalent nodal loading
    end
    elemENL(elnum,:)=elemENL(elnum,:)+ENL;              %Saving each member equivalent nodal loads to one matrix containing all member equivalent nodal loads
end

colENL=zeros(DOF*size(Node,1),1);
cellelemENL=cell(1,size(Member,1));
for elnum=1:size(elemENL,1)
    inode=Member(elnum,2);
    jnode=Member(elnum,3);
    index = [DOF*inode-((DOF-1):-1:0), DOF*jnode-((DOF-1):-1:0)];
    colelemENL=elemENL(elnum,2:end)';
    cellelemENL{elnum}=colelemENL;
    colENL(index)=colENL(index)+colelemENL;
end

nforce=zeros(DOF*size(Node,1),1);
for i=1:size(P,1)
    index=DOF*P(i,1)-((DOF-1):-1:0);    %Indexing of each nodal loads
    nforce(index)=P(i,2:DOF+1);         %Nodal load values for each DOF
end
force=colENL+nforce;                   %Net values of nodal loadings
%% Stiffness Assembly (K)
K = zeros(DOF*size(Node,1));    %Global stiffness, DOF*numberofnodes
for i=1:size(Member,1)          %Automate member-i local stiffness
    inode=Member(i,2);          %Identifying i node of each member
    jnode=Member(i,3);          %Identifying j node of each member
    xi = Node(inode,2);         %Identifying i node x position
    yi = Node(inode,3);         %Identifying i node y position
    xj = Node(jnode,2);         %Identifying j node x position
    yj = Node(jnode,3);         %Identifying j node y position
    index = [DOF*inode-((DOF-1):-1:0),DOF*jnode-((DOF-1):-1:0)]; %Indexing for each DOF (Example: Node 1 will have index (1 2), while node 4 will have index (7 8)
    L=sqrt((xj-xi)^2+(yj-yi)^2);%Member length
    c=(xj-xi)/L;
    s=(yj-yi)/L;
    
    MemberNumber=i;             %To show loop number (in this case equals member number)
    ElemType=Member(i,4);
    if isscalar(E)==1           %Taking account of differences in modulus of elasticity
        Eel=E;                  %If different modulus of elasticity exists for each member, iterate to each one
    else
        Eel=E(i,2);
    end
    if isscalar(I)==1           %Taking account of differences in sectional area
       Iel=I;                   %If different moment of inertia exists for each member, iterate to each one
    else
       Iel=I(i,2);
    end
    if isscalar(A)==1           %Taking account of differences in sectional area
       Ael=A;                   %If different sectional area exists for each member, iterate to each one
    else
       Ael=A(i,2);
    end
    [ki,T,k] = planeframe(Eel,Ael,Iel,L,c,s,ElemType); %Member's direct stiffness
    K(index,index) = K(index,index) + k;      %Save each looping values to each corresponding index
end
%% Boundary condition indexing
nDOFindex = (repmat(DOF*nDOF(:,1),1,DOF));      %Indexing of nDOF
for i=1:(DOF-1)
    nDOFindex(:,i)= nDOFindex(:,i)-(DOF-i);
end
nDOFindex = nDOFindex';
idnDOF  = logical(nDOF(:,2:end))';               %Turning nDOF value to logical (non zero to true, zero to false)
zeroDOF = nDOFindex(idnDOF);                     %Logical indexing (column space), taking only the restrained nDOF
freeDOF = setxor(zeroDOF,(1:DOF*size(Node,1)));  %Index of unrestrained DOF
%% Processing
reducedK = K(freeDOF,freeDOF);                  %Reducing K, obtaining only K components at non-zero DOF
reducedforce = force(freeDOF);                  %Reducing force, obtaining only forces at non-zero DOF
U=reducedK\reducedforce;                        %F=KU, U=inv(K)*F, thus obtaining displacements
Disp = zeros(DOF*size(Node,1),1);       
Disp(freeDOF) = U;                              %Displacements at indexed form
Displacement  = reshape(Disp,DOF,[])';          %Turning U column into u and v column
Displacements = [Node(:,1) Displacement];
DisplacementsTable=array2table(Displacements,'VariableNames',{'Node','u','v','theta'});
disp('<strong> ------------------ Displacements ------------------ </strong>')
disp('Note: u and v values in meter, theta values in radian')
disp(DisplacementsTable)

reducedK2 = K(zeroDOF,freeDOF);                 %Reducing K, taking restrained rows and unrestrained column
support = reducedK2*U-force(zeroDOF);           %Calculating forces at zeroDOFs, F=KU, obtaining support reactions
Reac = zeros(DOF*size(Node,1),1);
Reac(zeroDOF) = support;                        %Reactions at indexed form
Reaction = reshape(Reac,DOF,[])';               %Turning support column into FX and FY column
Reactions=[Node(:,1) Reaction];  
IndicesReaction=setxor(Reactions(:,1),nDOF(:,1));    
Reactions(IndicesReaction,:) = [];
ReactionsTable=array2table(Reactions,'VariableNames',{'Node','FX','FY','M'});
disp('<strong> ------------------ Support Reactions ------------------ </strong>')
disp('Note: FX and FY in kN, M in kNm')
disp(ReactionsTable)

Internal = zeros(size(Member,1),1+2*DOF);
Internal(:,1) = 1:size(Member,1);
for elnum=1:size(Member,1)
    inode=Member(elnum,2);       %Identifying i node of each member
    jnode=Member(elnum,3);       %Identifying j node of each member
    xi = Node(inode,2);          %Identifying i node x position
    yi = Node(inode,3);          %Identifying i node y position
    xj = Node(jnode,2);          %Identifying j node x position
    yj = Node(jnode,3);          %Identifying j node y position
    L=sqrt((xj-xi)^2+(yj-yi)^2);%Member length
    c=(xj-xi)/L;
    s=(yj-yi)/L;
    
    MemberNumber=elnum;         %To show loop number (in this case equals member number)
    if isscalar(E)==1           %Taking account of differences in modulus of elasticity
        Eel=E;                  %If different modulus of elasticity exists for each member, iterate to each one
    else
        Eel=E(elnum,2);
    end
    if isscalar(I)==1           %Taking account of differences in moment of inertia
       Iel=I;                   %If different moment of inertia exists for each member, iterate to each one
    else
       Iel=I(elnum,2);
    end
    if isscalar(A)==1           %Taking account of differences in sectional area
       Ael=A;                   %If different cross sectional area exists for each member, iterate to each one
    else
       Ael=A(elnum,2);
    end
    
    index = [DOF*inode-((DOF-1):-1:0),DOF*jnode-((DOF-1):-1:0)];
    ElemType=Member(elnum,4);
    [ki,T,k] = planeframe(Eel,Ael,Iel,L,c,s,ElemType);
    elemInternal = ki*T*Disp(index)-T'\cellelemENL{elnum};
    Internal(elnum,2:7)=(elemInternal'+Internal(elnum,2:7));
    Axial=Internal(:,[1 2 5]);
    Shear=Internal(:,[1 3 6]);
    BM=Internal(:,[1 4 7]);
end

AxialTable=array2table(Axial,'VariableNames',{'Member','iaxial','jaxial'});
disp('<strong> ------------------ Axial Forces ------------------ </strong>')
disp('Note: Values in kN')
disp(AxialTable)

ShearTable=array2table(Shear,'VariableNames',{'Member','ishear','jshear'});
disp('<strong> ------------------ Shear Forces ------------------ </strong>')
disp('Note: Values in kN')
disp(ShearTable)

BMTable=array2table(BM,'VariableNames',{'Member','imoment','jmoment'});
disp('<strong> ------------------ Bending Moments ------------------ </strong>')
disp('Note: Values in kNm')
disp(BMTable)
%% Figure Plots
figplot=figure();
figplot.WindowState='maximized';
for Thp=1:size(Member,1)
    if Member(Thp,4) == 2
    xundef=[Node(Member(Thp,2),2) Node(Member(Thp,3),2)];
    yundef=[Node(Member(Thp,2),3) Node(Member(Thp,3),3)];
    subplot(3,3,[1,2,4,5,7,8]);
    Pundef=plot(xundef,yundef,'.-b','linewidth',0.5);xlim([min(Node(:,2))-1 max(Node(:,2))+1]),ylim([min(Node(:,3))-1 max(Node(:,3))+1]); xlabel('X Coordinates [m]'), ylabel('Y Coordinates [m]'), xticks(0:1:max(Node(:,2))), grid('on')
    Pundef.MarkerSize=10;
    hold on
    xel=(xundef(1,1)+xundef(1,2))/2;yel=(yundef(1,1)+yundef(1,2))/2;
    text(xel, yel, int2str(Thp), 'fontsize', 10);
    else
    xundef=[Node(Member(Thp,2),2) Node(Member(Thp,3),2)];
    yundef=[Node(Member(Thp,2),3) Node(Member(Thp,3),3)];
    subplot(3,3,[1,2,4,5,7,8]);
    Pundef=plot(xundef,yundef,'.-b','linewidth',1.5);xlim([min(Node(:,2))-1 max(Node(:,2))+1]),ylim([min(Node(:,3))-1 max(Node(:,3))+1]); xlabel('X Coordinates [m]'), ylabel('Y Coordinates [m]'), xticks(0:1:max(Node(:,2))), grid('on')
    Pundef.MarkerSize=10;
    hold on
    xel=(xundef(1,1)+xundef(1,2))/2;yel=(yundef(1,1)+yundef(1,2))/2;
    text(xel, yel, int2str(Thp), 'fontsize', 10);
    end
NNode = [Node(:,1) Displacement(:,1:2)*FigScale+Node(:,2:3)];
    xdef=[NNode(Member(Thp,2),2) NNode(Member(Thp,3),2)];
    ydef=[NNode(Member(Thp,2),3) NNode(Member(Thp,3),3)];
subplot(3,3,[1,2,4,5,7,8]);
plot(xdef,ydef,'-.r');title(['Chord Movements (Scaled ', num2str(FigScale), 'x)'])
hold on
legend('Undeformed','Deformed','location','southeast')
end
for i=1:size(Member,1)
    inode=Member(i,2);
    jnode=Member(i,3);
    ix=Node(inode,2);
    iy=Node(inode,3);
    jx=Node(jnode,2);
    jy=Node(jnode,3);
    l=sqrt((jx-ix)^2+(jy-iy)^2);
    c=(jx-ix)/l;
    s=(jy-iy)/l;
    vi=Shear(i,2);
    vj=-Shear(i,3);
    dxvi=ix-vi*s*ShearScale;
    dyvi=iy+vi*c*ShearScale;
    dxvj=jx-vj*s*ShearScale;
    dyvj=jy+vj*c*ShearScale;
    ShearDiagramx=[dxvi dxvj];
    ShearDiagramy=[dyvi dyvj];
    subplot(3,3,3)
    ppplot=plot(ShearDiagramx,ShearDiagramy,'-g');
    hold on;
    undefx=[ix jx];
    undefy=[iy jy];
    if Member(i,4)==1
    Pundef=plot(undefx,undefy,'.-b','linewidth',1.5);xlim([min(Node(:,2))-1 max(Node(:,2))+1]),ylim([min(Node(:,3))-1 max(Node(:,3))+1]),grid('on')
    Pundef.MarkerSize= 15;
    else
    Pundef=plot(undefx,undefy,'.-b','linewidth',0.5);xlim([min(Node(:,2))-1 max(Node(:,2))+1]),ylim([min(Node(:,3))-1 max(Node(:,3))+1]),grid('on')
    Pundef.MarkerSize= 15;
    end
    hold on;
    txi=[dxvi ix];
    tyi=[dyvi iy];
    plot(txi,tyi,'-g');
    hold on;
    txj=[dxvj jx];
    tyj=[dyvj jy];
    plot(txj,tyj,'-g');
    title(['Shear Force Diagram (Scaled ', num2str(ShearScale), 'x)'])
    hold on;
end
for i=1:size(Member,1)
    inode=Member(i,2);
    jnode=Member(i,3);
    ix=Node(inode,2);
    iy=Node(inode,3);
    jx=Node(jnode,2);
    jy=Node(jnode,3);
    l=sqrt((jx-ix)^2+(jy-iy)^2);
    c=(jx-ix)/l;
    s=(jy-iy)/l;
    mi=BM(i,2);
    mj=-BM(i,3);
    dxmi=ix-mi*s*BendingMomentScale;
    dymi=iy+mi*c*BendingMomentScale;
    dxmj=jx-mj*s*BendingMomentScale;
    dymj=jy+mj*c*BendingMomentScale;
    MomentDiagramx=[dxmi dxmj];
    MomentDiagramy=[dymi dymj];
    subplot(3,3,6);
    plot(MomentDiagramx,MomentDiagramy,'-r');
    hold on;
    undefx=[ix jx];
    undefy=[iy jy];
    if Member(i,4)==1
    Pundef=plot(undefx,undefy,'.-b','linewidth',1.5);xlim([min(Node(:,2))-1 max(Node(:,2))+1]),ylim([min(Node(:,3))-1 max(Node(:,3))+1]),grid('on')
    Pundef.MarkerSize= 15;
    else
    Pundef=plot(undefx,undefy,'.-b','linewidth',0.5);xlim([min(Node(:,2))-1 max(Node(:,2))+1]),ylim([min(Node(:,3))-1 max(Node(:,3))+1]),grid('on')
    Pundef.MarkerSize= 15;
    end
    hold on;
    txi=[dxmi ix];
    tyi=[dymi iy];
    plot(txi,tyi,'-r');
    hold on;
    txj=[dxmj jx];
    tyj=[dymj jy];
    plot(txj,tyj,'-r');
    hold on;
    title(['Bending Moment Diagram (Scaled ', num2str(BendingMomentScale), 'x)'])
end
for i=1:size(Member,1)
    inode=Member(i,2);
    jnode=Member(i,3);
    ix=Node(inode,2);
    iy=Node(inode,3);
    jx=Node(jnode,2);
    jy=Node(jnode,3);
    l=sqrt((jx-ix)^2+(jy-iy)^2);
    c=(jx-ix)/l;
    s=(jy-iy)/l;
    ai=-Axial(i,2);
    aj=Axial(i,3);
    dxai=ix-ai*s*AxialScale;
    dyai=iy+ai*c*AxialScale;
    dxaj=jx-aj*s*AxialScale;
    dyaj=jy+aj*c*AxialScale;
    AxialDiagramx=[dxai dxaj];
    AxialDiagramy=[dyai dyaj];
    subplot(3,3,9);
    plot(AxialDiagramx,AxialDiagramy,'-c');
    hold on;
    undefx=[ix jx];
    undefy=[iy jy];
    if Member(i,4)==1
    Pundef=plot(undefx,undefy,'.-b','linewidth',1.5);xlim([min(Node(:,2))-1 max(Node(:,2))+1]),ylim([min(Node(:,3))-1 max(Node(:,3))+1]),grid('on')
    Pundef.MarkerSize= 15;
    else
    Pundef=plot(undefx,undefy,'.-b','linewidth',0.5);xlim([min(Node(:,2))-1 max(Node(:,2))+1]),ylim([min(Node(:,3))-1 max(Node(:,3))+1]),grid('on')
    Pundef.MarkerSize= 15;
    end
    hold on;
    txi=[dxai ix];
    tyi=[dyai iy];
    plot(txi,tyi,'-c');
    hold on;
    txj=[dxaj jx];
    tyj=[dyaj jy];
    plot(txj,tyj,'-c');
    hold on;
    title(['Axial Force Diagram (Scaled ', num2str(AxialScale), 'x)'])
end
%% Functions
function [ki,T,k] = planeframe(Eel,Ael,Iel,L,c,s,ElemType)
if ElemType == 2
ki  = [Ael*Eel/L 0 0 -Ael*Eel/L 0 0
      0 0 0	0 0 0
      0 0 0 0 0 0
      -Ael*Eel/L 0 0 Ael*Eel/L 0 0
      0 0 0 0 0 0
      0 0 0 0 0 0];
elseif ElemType == 1
ki = [Ael*Eel/L 0 0 -Ael*Eel/L 0 0
      0 12*Eel*Iel/L^3 6*Eel*Iel/L^2 0 -12*Eel*Iel/L^3 6*Eel*Iel/L^2
      0 6*Eel*Iel/L^2 4*Eel*Iel/L 0 -6*Eel*Iel/L^2 2*Eel*Iel/L
      -Ael*Eel/L 0 0 Ael*Eel/L 0 0
      0 -12*Eel*Iel/L^3 -6*Eel*Iel/L^2 0 12*Eel*Iel/L^3 -6*Eel*Iel/L^2
      0 6*Eel*Iel/L^2 2*Eel*Iel/L 0 -6*Eel*Iel/L^2 4*Eel*Iel/L];
end
  
T = [c s 0 0 0 0
     -s c 0 0 0 0
     0 0 1 0 0 0
     0 0 0 c s 0
     0 0 0 -s c 0
     0 0 0 0 0 1];
 
k = T'*ki*T;
end
function [enli,T,enl]=interload(w1t,w2t,w1a,w2a,l,c,s)
enli=[-l*(2*w1a+w2a)/6 -(3*w2t+7*w1t)*l/20 -(3*w1t+2*w2t)*l^2/60 -l*(w1a+2*w2a)/6 -(3*w1t+7*w2t)*l/20 (2*w1t+3*w2t)*l^2/60];

T = [c s 0 0 0 0
     -s c 0 0 0 0
     0 0 1 0 0 0
     0 0 0 c s 0
     0 0 0 -s c 0
     0 0 0 0 0 1];
 
enl = T'*enli';
end