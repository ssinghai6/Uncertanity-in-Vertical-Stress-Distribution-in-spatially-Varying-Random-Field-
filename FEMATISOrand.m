function [matmtrx, Elas1, elas]=FEMATISOrand(iopt,poisson,cl,cov,emean)

%------------------------------------------------------------------------
%  Purpose:
%     determine the constitutive equation for isotropic material
%
%  Synopsis:
%     [matmtrx]=fematiso(iopt,elastic,poisson) 
%
%  Variable Description:
%     elastic - elastic modulus
%     poisson - Poisson's ratio   
%     iopt=1 - plane stress analysis
%     iopt=2 - plane strain analysis
%     iopt=3 - axisymmetric analysis
%     iopt=4 - three dimensional analysis
%      n- no. of elements
%     Emod- Elastic modulus for generating random field
%     a(:,:,1) or a(:,:,2) will form as multi-dimensional array
%------------------------------------------------------------------------
n=3960;                  % number of elements
nnel=8;                  % number of nodes per element
ndof=3;                  % number of dofs per node
nnode=4807;              % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
%poisson=0.3;             % Poisson's ratio

nglx=2; ngly=2;          % 2x2 Gauss-Legendre quadrature
nglxy=nglx*ngly;         % number of sampling points per element
%cl=1.5;                    % Markov correlation length
%emean=50;                % Mean Elastic Modulus
%cov=30;                  % coeff. of variation in (%)
[Elas1,elas]=randomf3dcmd(cl,cov,emean);

for i=1:n
 if iopt==1        % plane stress
   matmtrx(:,:,i)= Elas1(i)/(1-poisson*poisson)* ...
   [1  poisson 0; ...
   poisson  1  0; ...
   0  0  (1-poisson)/2];

 elseif iopt==2     % plane strain
   matmtrx(:,:,i)= Elas1(i)/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson 0; 
   poisson  (1-poisson)  0;
   0  0  (1-2*poisson)/2];

 elseif iopt==3     % axisymmetry
   matmtrx(:,:,i)= Elas1(i)/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson  0; 
   poisson  (1-poisson)   poisson  0;
   poisson  poisson  (1-poisson)   0;
   0    0    0   (1-2*poisson)/2];
 
 else     % three-dimension
   matmtrx(:,:,i)= Elas1(i)/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];

 end
end
