cl=[2.5];
cov=[30]; 
emean=[50];
nnode=4807;
nit=200;
j=length(cl);
k=length(cov);
l=length(emean);
sigmay = zeros(4807, nit, l*j*k);
dispy1 = zeros(4807, nit, l*j*k);
Elas1 = zeros(3960, nit, l*j*k);
elas = zeros(4807, nit, l*j*k);
count = 0;
for a = cl;
for b = cov
    count = count + 1
for i=1:nit
    [sigmay(:,i, count),dispy1(:,i, count),Elas1(:,i, count),elas(:,i, count),norm(:,i, count)]= FinalrandomFEM3D(cl,b,emean);
end
end
end
for i = 1:count;
a1= 'dispy1';
b1= 'Elas1';
c1= 'elas';
d1= 'norm';
e1= 'sigmay';
a2 = num2str(i);
a3 = '.xlsx';
val1 = strcat(a1, a2, a3);
val2 = strcat(b1, a2, a3);
val3 = strcat(c1, a2, a3);
val4 = strcat(d1, a2, a3);
val5 = strcat(e1, a2, a3);

xlswrite(val1, dispy1(:, :, i));
xlswrite(val2, Elas1(:, :, i));
xlswrite(val3, elas(:, :, i));
xlswrite(val4, norm(:, :, i));
xlswrite(val5, sigmay(:, :, i));
end