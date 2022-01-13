   


%----------------------------------------------------------------------------
%                                                              
%   3D Analysis
%
% Variable descriptions                                                      
%   k = element matrix                                            
%   f = element vector
%   kk = system matrix                                            
%   ff = system vector                                                
%   disp = system nodal displacement vector
%   eldisp = element nodal displacement vector
%   stress = matrix containing stresses
%   strain = matrix containing strains
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element    
%   point2 = matrix containing sampling points
%   weight2 = matrix containing weighting coefficients
%   bcdof = a vector containing dofs associated with boundary conditions    
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------
function[sigmay,dispy1,Elas1,elas,norm]=FinalrandomFEM3D(cl,cov,emean)

nel=3960;                % number of elements
nnel=8;                  % number of nodes per element
ndof=3;                  % number of dofs per node
nnode=4807;              % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
%emodule=100;              % elastic modulus
poisson=0.3;             % Poisson's ratio
nglx=2; ngly=2; nglz =2;         % 3x3 Gauss-Legendre quadrature
nglxy=nglx*ngly*nglz;         % number of sampling points per element
               
%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------
load nodes.txt
load elements.txt
load bound.txt
A=nodes;
B=elements;

gcoord = zeros(nnode,2);

 gcoord =A;
 node= B;
  %-------------------------------------
%  input data for boundary conditions
%-------------------------------------
bcdof=bound';
 
 
bcval= zeros(1,length(bcdof)) ;      

       

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
disp1=zeros(sdof,1);     % system displacement vector
eldisp=zeros(edof,1);   % element displacement vector
stress=zeros(nglxy,6);  % matrix containing stress components
strain=zeros(nglxy,6);  % matrix containing strain components
index=zeros(edof,1);    % index vector
kinmtx3=zeros(6,edof);   % kinematic matrix
matmtx=zeros(6,6);      % constitutive matrix
finalstress = zeros(nnel,6,nel);
finalstress1 = zeros(nnel,6,nel);
extp = zeros(nnel,nnel);
%----------------------------
%  force vector
%----------------------------
%conversion of udl into nodal loads
%load is acting in y-direction
%node nubers on which udl is acting in a sequence

%position of nodal load into system force vector
ff(1934) =-6.25;
ff(2006) =-12.5;
ff(2069) =-12.5;
ff(2207) =-12.5;
ff(2405) =-12.5;
ff(2522) =-12.5;
ff(2576) =-25;
ff(2630) =-12.5;
ff(2648) =-25;
ff(2795) =-25;
ff(2963) =-12.5;
ff(2999) =-25;
ff(3128) =-6.25;
ff(3176) =-12.5;
ff(3230) =-25;
ff(3269) =-12.5;
ff(3341) =-12.5;
ff(3431) =-12.5;
ff(3542) =-25;
ff(3587) =-12.5;
ff(3725) =-12.5;
ff(3830) =-12.5;
ff(3899) =-25;
ff(4154) =-12.5;
ff(4205) =-12.5;
ff(4334) =-25;
ff(4472) =-12.5;
ff(4685) =-6.25;
ff(4775) =-25;
ff(4919) =-12.5;
ff(5249) =-12.5;
ff(5342) =-12.5;
ff(5891) =-6.25;
qo = sum(ff)/5; %load applied
%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

[point3,weight3]=FEGLQD3(nglx,ngly,nglz);       % sampling points & weights
[matmtx, Elas1, elas]=FEMATISOrand(4,poisson,cl,cov,emean);        % compute constitutive matrix

for iel=1:nel           % loop for the total number of elements

for i=1:nnel
nd(i)=node(iel,i);         % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);   % extract y value of the node
zcoord(i)=gcoord(nd(i),3);   % extract z value of the node
end

k=zeros(edof,edof);         % initialization of element matrix to zero

%--------------------------------
%  numerical integration
%--------------------------------
for intx = 1:nglx
    x = point3(intx,1);
    wtx = weight3(intx,1);
for inty = 1:ngly
    y = point3(inty,1);
    wty = weight3(inty,1);
for intz = 1:nglz
    z = point3(intz,1);
    wtz = weight3(intz,1);


[shape,dhdr,dhds,dhdt]=FEISOS8(x,y,z);    % compute shape functions and
                                   % derivatives at sampling point

jacob3=FEJACOB3(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord);  % compute Jacobian

detjacob=det(jacob3);                 % determinant of Jacobian
invjacob=inv(jacob3);                 % inverse of Jacobian matrix

[dhdx,dhdy,dhdz]=FEDERIV3(nnel,dhdr,dhds,dhdt,invjacob); % derivatives w.r.t.
                                               % physical coordinate

kinmtx3=FEKINE3D(nnel,dhdx,dhdy,dhdz);          % compute kinematic matrix

%------------------------------
%  compute element matrix
%------------------------------

k=k+kinmtx3'*matmtx(:,:,iel)*kinmtx3*wtx*wty*wtz*detjacob;    % element matrix

end
end
end
                                   % end of numerical integration loop

index=FEELDOF(nd,nnel,ndof);% extract system dofs associated with element

kk=FEASMBL1(kk,k,index);  % assemble element matrices

end
Atemp = kk;
ftemp = ff;
index = length(bound);
for i = sdof:-1:1
    if i == bound(index)
        index = index-1;
        Atemp(:, i) = [];
        Atemp(i, :) = [];
        ftemp(i) = [];
    end
end
       
%----------------------------
%  solve the matrix equation
%----------------------------

disptemp = Atemp\ftemp;
disp1 = zeros(sdof,1);
temp = 0;
index = 1;
sumTemp = 0;
for i = 1:sdof
    if i == bound(index)
        index = index+1;
        sumTemp = sumTemp+1;
        disp1(i) = 0;
    else
        temp = temp+1;
        disp1(i) = disptemp(temp);
    end
end  


for i=1:nnode
    dispy1(i)=disp1((3*i)-1);
end
%stress calculation
for ielp=1:nel           % loop for the total number of elements

for i=1:nnel
nd(i)=node(ielp,i);        % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2); % extract y value of the node
zcoord(i)=gcoord(nd(i),3); % extract z value
end

%--------------------------------
%  numerical integration
%--------------------------------

intp=0;

for intx = 1:nglx
    x = point3(intx,1);
    wtx = weight3(intx,1);
for inty = 1:ngly
    y = point3(inty,1);
    wty = weight3(inty,1);
for intz = 1:nglz
    z = point3(intz,1);
    wtz = weight3(intz,1);
 intp = intp+1;


[shape,dhdr,dhds,dhdt]=FEISOS8(x,y,z);   % compute shape functions and
                                   % derivatives at sampling point

jacob3=FEJACOB3(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord);  % compute Jacobian

detjacob=det(jacob3);                 % determinant of Jacobian
invjacob=inv(jacob3);                 % inverse of Jacobian matrix

[dhdx,dhdy,dhdz]=FEDERIV3(nnel,dhdr,dhds,dhdt,invjacob); % derivatives w.r.t.
                                               % physical coordinate

kinmtx3=FEKINE3D(nnel,dhdx,dhdy,dhdz);     % compute kinematic matrix
index=FEELDOF(nd,nnel,ndof);% extract system dofs for the element

%-------------------------------------------------------
%  extract element displacement vector
%-------------------------------------------------------

for i=1:edof
eldisp(i)=disp1(index(i));
end        
estrain=kinmtx3*eldisp;             % compute strains
estress=matmtx(:,:,ielp)*estrain;             % compute stresses

for i=1:6
strain(intp,i)=estrain(i);          % store for each element
stress(intp,i)=estress(i);      % store for each element          
end
A = (5+(3*sqrt(3)))/4;
B =-(sqrt(3)+1)/4;
C = (sqrt(3)-1)/4;
D = (5-(3*sqrt(3)))/4;
extp(1,:)= [A B C B B C D C];
extp(2,2:8)= [A B C C B C D];
extp(3,3:8)=[A B D C B C];
extp(4,4:8)=[A C D C B];
extp(5,5:8)=[A B C B];
extp(6,6:8)=[A B C];
extp(7,7:8)=[A B];
extp(8,8)=[A];
for i = 1:8
    for j = 1:8
  extp(j,i) = extp(i,j);
    end
  %extp = ones(8,8);
location=[ielp,intx,inty,intz] ;        % print location for stress
stress(intp,:);                    % print stress values
if intx==2&&inty==2
   num2=stress;
   finalstress1(:,:,ielp)=num2;
    finalstress(:,:,ielp)=extp*num2;
    
   
   
end
end
                                 


end

for ind = 1:nnode
    [r c] = find(node==ind) ;
    share = length(r) ;
   
    avg_stress2(ind,1)=ind;  
    if share ==1
        avg_stress2(ind,2)=finalstress(c,2,r);
        disp (finalstress(c,2,r));
   
    else if share ==2
     avg_stress2(ind,2)=(finalstress(c(1),2,r(1))+finalstress(c(2),2,r(2)))/2;
   else if share ==3
       avg_stress2(ind,2)=(finalstress(c(1),2,r(1))+finalstress(c(2),2,r(2))+finalstress(c(3),2,r(3)))/3;
        else if share ==4
        avg_stress2(ind,2)= (finalstress(c(1),2,r(1))+finalstress(c(2),2,r(2))+finalstress(c(3),2,r(3))+finalstress(c(4),2,r(4)))/4;
        else if share ==5
        avg_stress2(ind,2)= (finalstress(c(1),2,r(1))+finalstress(c(2),2,r(2))+finalstress(c(3),2,r(3))+finalstress(c(4),2,r(4))+finalstress(c(5),2,r(5)))/5;
         else if share ==6
        avg_stress2(ind,2)=  (finalstress(c(1),2,r(1))+finalstress(c(2),2,r(2))+finalstress(c(3),2,r(3))+finalstress(c(4),2,r(4))+finalstress(c(5),2,r(5))+finalstress(c(6),2,r(6)))/6;
         else if share ==7
        avg_stress2(ind,2)= (finalstress(c(1),2,r(1))+finalstress(c(2),2,r(2))+finalstress(c(3),2,r(3))+finalstress(c(4),2,r(4))+finalstress(c(5),2,r(5))+finalstress(c(6),2,r(6))+finalstress(c(7),2,r(7)))/7;
         else if share ==8
        avg_stress2(ind,2)= (finalstress(c(1),2,r(1))+finalstress(c(2),2,r(2))+finalstress(c(3),2,r(3))+finalstress(c(4),2,r(4))+finalstress(c(5),2,r(5))+finalstress(c(6),2,r(6))+finalstress(c(7),2,r(7))+finalstress(c(8),2,r(8)))/8;
            end
             end
             end
             end
            end
            end
        end
   end
end
end
sigmay=avg_stress2(:,2);
norm=sigmay(:,1)/qo;

end
end
end