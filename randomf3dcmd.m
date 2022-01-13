%**********************************************************
%Generation of random 3D Random field using Correlation Matrix Decomposition
%single random paramete
function[Elas1,elas]=randomf3dcmd(cl,cov,emean)
% Elas for elements and elas is for nodes
load cc.txt
load elements.txt
%element txt is the node connectivity matrix
%cc is txt file containing nodes(sequence,x,y)
n=size(cc(:,1));
%cl=5;% cl is the correlation length
%cov=10; %coefficient of variation( for 50% cov=50)
%emean=50;%(mean value of elastic modulus in Mpa)
for i=1:n
for j=1:n
d11(i,j)=sqrt((cc(j,4)-cc(i,4))^2+(cc(j,3)-cc(i,3))^2+((cc(j,2)-cc(i,2))^2));
%d11 is the distance between two nodes (n*n)
e1(i,j)=exp(-2*d11(i,j)/cl);

%markov correlation function where d11 is distance
%in denominator we can vary degree of fluctuation
%decomposition of e1(markov) into lower and upper triangle usings
%cholskey's decomposition(for this matrix should be a square matrix
end
end
[L,U] = lu(e1);

n=length(cc(:,1));
ne=length(elements(:,1));
%cov=50;
%emean=50;
mu = 0;
sigma= 1;
b11= normrnd(mu,sigma,n,1);
%generation of random no. b1(n*1) which are normally distributed with
%sigma=1 and mu=0

mm = emean;
vv = (mm*cov)/100;
sigmaf = sqrt(log((cov/100)^2+1));
%muf = log((mm^2)/sqrt(vv+mm^2));
muf=log(mm)-(0.5*((sigmaf)^2));
%mm is mean and vv is sigma for desired set of generated e values
%cov=mu/sigma(cov varies from 5% to 50%)

Z=L*b11;

%z(n*1)=L(n*n)*b11(n*1)(L formed by cholkey's decomposition and b11 are randomly generate no.) 
for k=1:n
    elas(k,1)=exp(muf +(sigmaf*Z(k,1)));
    %elas1 (n*1)- elastic mod with respect to nodes.    
    
end
%Elastic modulus calculation for elements ie, taking average of emod at nodes
c=zeros(n,5);
c(:,1:4)=cc(:,1:4);
c(:,5)=elas(:,1);
%load c.txt
%load join.txt
a = zeros(1,ne);
c(1:n,5)=elas(1:n,1);
for index = 1:ne
    sum1 = 0;
    for r = elements(index,1:8)
        sum1=sum1+c(r, 5);
    end
    a(index)=sum1/8;
end
Elas1=a';
end