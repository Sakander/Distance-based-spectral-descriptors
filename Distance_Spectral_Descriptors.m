close all
clc
format short g
roundn = @(t,n) round(t*10^n)./10^n;

n=14;            % Order of the graph

B=[[0, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 5, 6, 5], [1, 0, 2, 1, 3, 2, 3, 4, 5, 4, 5, 6, 7, 6], [1, 2, 0, 3, 1, 2, 1, 2, 3, 2, 3, 4, 5, 4], [2, 1, 3, 0, 2, 1, 4, 5, 4, 3, 6, 7, 6, 5], [2, 3, 1, 2, 0, 1, 2, 3, 2, 1, 4, 5, 4, 3], [3, 2, 2, 1, 1, 0, 3, 4, 3, 2, 5, 6, 5, 4], [2, 3, 1, 4, 2, 3, 0, 1, 2, 3, 2, 3, 4, 3], [3, 4, 2, 5, 3, 4, 1, 0, 1, 2, 1, 2, 3, 2], [4, 5, 3, 4, 2, 3, 2, 1, 0, 1, 2, 3, 2, 1], [3, 4, 2, 3, 1, 2, 3, 2, 1, 0, 3, 4, 3, 2], [4, 5, 3, 6, 4, 5, 2, 1, 2, 3, 0, 1, 2, 3], [5, 6, 4, 7, 5, 6, 3, 2, 3, 4, 1, 0, 1, 2], [6, 7, 5, 6, 4, 5, 4, 3, 2, 3, 2, 1, 0, 1], [5, 6, 4, 5, 3, 4, 3, 2, 1, 2, 3, 2, 1, 0]];



D= reshape(B,[n,n]);
if D==D'
      disp('Matrix is Symmetric')
else
      disp('Not Symmetric')
end

% The distance matrix

q=eig(D);
x = roundn(q,4);
rho=max(x);
E=sum(abs(eig(D)));
EE=sum(exp(eig(D)));

% The distance Laplacian matrix
Tr=sum(D)';
t=sum(Tr)/n;
sig=sum(Tr)/2;
Diag=diag(Tr);
DL=Diag-D;
q1=eig(DL);
x1 = roundn(q1,4);
rho1=max(x1);
E1=sum(abs(eig(DL)-t));
EE1=roundn(sum(exp(eig(DL)-sig)),4);

% The distance signless Laplacian matrix

DQ=Diag+D;
q2=eig(DQ);
x2 = roundn(q2,4);
rho2=max(x2);
E2=sum(abs(eig(DQ)-t));
EE2=roundn(sum(exp(eig(DQ)-sig)),4);

% The Schultz matrix
diam=max(max(D));
A=[];
for r=2:diam
    D(D==r)=0;
    A=D;
end
k=sum(A)';
D1= reshape(B,[n,n]);
S=A+D1;
q3=eig(S);
x3 = roundn(q3,4);
rho3=max(x3);
E3=sum(abs(eig(S)));
EE3=sum(exp(eig(S)));

% The Harary matrix
l=size(D1,1);
H=[];
for i=1:l
    for j=1:l
        if i~=j
            D1(i,j)=1/D1(i,j);
            H=D1;
        end
    end
end
q4=eig(H);
x4 = roundn(q4,4);
rho4=max(x4);
E4=sum(abs(eig(H)));
EE4=sum(exp(eig(H)));

% The degree-distance matrix
D2= reshape(B,[n,n]);
DD=[];
for i=1:l
    for j=1:l
        if i~=j
            D2(i,j)=(k(i)+k(j))/D2(i,j);
            DD=D2;
        end
    end
end
q5=eig(DD);
x5 = roundn(q5,4);
rho5=max(x5);
E5=sum(abs(eig(DD)));
EE5=sum(exp(eig(DD)));

% The Gutman matrix
D3= reshape(B,[n,n]);
Gut=[];
for i=1:l
    for j=1:l
        if i~=j
            D3(i,j)=(k(i)*k(j))/D3(i,j);
            Gut=D3;
        end
    end
end
q6=eig(Gut);
x6 = roundn(q6,4);
rho6=max(x6);
E6=sum(abs(eig(Gut)));
EE6=sum(exp(eig(Gut)));


% The Szeged matrix
D4= reshape(B,[n,n]);
l=size(D4,1);
t1=[];
s1=[];
e=[];
for i=1:l
    for j=i+1:l
        if D4(i,j)==1
            e=[e;[i,j]];
        end
    end
end

ress = zeros(l,l);
rest = zeros(l,l);
resk = zeros(l,l);
for i=1:size(e,1)
    ss = [];
    tt = [];
    kk = [];
    for ii=1:l
        if D4(e(i,1),ii)<D4(e(i,2),ii)
           ss=[ss ii];
        end
        if D4(e(i,1),ii)>D4(e(i,2),ii)
           tt=[tt ii];
         end
        if D4(e(i,1),ii)==D4(e(i,2),ii)
           kk=[kk ii];
        end
    end
    ress(i, 1:length(ss)) = ss;
    rest(i, 1:length(tt)) = tt;
    resk(i, 1:length(kk)) = kk;
end
res1=[e ress];
res2=[e rest];
res3=[e resk];
[sm, ~] = size(s1);
U=[];
for i=1:sm
    U = [U; length(s1(i,:))];
end
Us = (ress~=0);
U=sum(Us');
U=U';

[tm, ~] = size(t1);
V=[];
for i=1:tm
    V = [V; length(t1(i,:))];
end
Vt = (rest~=0);
V=sum(Vt');
V= V';
N=[U V];
Sz=[];
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=(N(i)*N(j));
            Sz=A;
        end
    end
end
q7=eig(Sz);
x7 = roundn(q7,4);
rho7=max(x7);
E7=sum(abs(eig(Sz)));
EE7=sum(exp(eig(Sz)))


% The PI matrix
diam=max(max(D4));
A1=[];
for r=2:diam
    D4(D4==r)=0;
    A1=D4;
end
PI=[];
for i=1:l
    for j=1:l
        if A1(i,j)==1
            A1(i,j)=(N(i)+N(j));
            PI=A1;
        end
    end
end
q8=eig(PI);
x8 = roundn(q8,4);
rho8=max(x8);
E8=sum(abs(eig(PI)));
EE8=sum(exp(eig(PI)));

% The ABC2 matrix
D5= reshape(B,[n,n]);
diam=max(max(D5));
A2=[];
for r=2:diam
    D5(D5==r)=0;
    A2=D5;
end
ABC2=[];
for i=1:l
    for j=1:l
        if A2(i,j)==1
            A2(i,j)=sqrt((N(i)+N(j)-2)/(N(i)*N(j)));
            ABC2=A2;
        end
    end
end
q9=eig(ABC2);
x9 = roundn(q9,4);
rho9=max(x9);
E9=sum(abs(eig(ABC2)));
EE9=sum(exp(eig(ABC2)));


% The GA2 matrix
D6= reshape(B,[n,n]);
diam=max(max(D6));
A3=[];
for r=2:diam
    D6(D6==r)=0;
    A3=D6;
end
GA2=[];
for i=1:l
    for j=1:l
        if A3(i,j)==1
            A3(i,j)=(2*sqrt(N(i)*N(j)))/((N(i)+N(j)));
            GA2=A3;
        end
    end
end
q10=eig(GA2);
x10 = roundn(q10,4);
rho10=max(x10);
E10=sum(abs(eig(GA2)));
EE10=sum(exp(eig(GA2)));


disp('Spectrum-related descriptors:')
fprintf('The D-SR is %4.4f\n',rho');
fprintf('The D-E is %4.4f\n',E');
fprintf('The D-EE is %4.4f\n',EE');
fprintf('The DL-SR is %4.4f\n',rho1');
fprintf('The DL-E is %4.4f\n',E1');
fprintf('The DL-EE is %4.4f\n',EE1');
fprintf('The DQ-SR is %4.4f\n',rho2');
fprintf('The DQ-E is %4.4f\n',E2');
fprintf('The DQ-EE is %4.4f\n',EE2');
fprintf('The S-SR is %4.4f\n',rho3');
fprintf('The S-E is %4.4f\n',E3');
fprintf('The S-EE is %4.4f\n',EE3');
fprintf('The H-SR is %4.4f\n',rho4');
fprintf('The H-E is %4.4f\n',E4');
fprintf('The H-EE is %4.4f\n',EE4');
fprintf('The DD-SR is %4.4f\n',rho5');
fprintf('The DD-E is %4.4f\n',E5');
fprintf('The DD-EE is %4.4f\n',EE5');
fprintf('The Gut-SR is %4.4f\n',rho6');
fprintf('The Gut-E is %4.4f\n',E6');
fprintf('The Gut-EE is %4.4f\n',EE6');
fprintf('The Sz-SR is %4.4f\n',rho7');
fprintf('The Sz-E is %4.4f\n',E7');
fprintf('The Sz-EE is %4.4f\n',EE7');
fprintf('The PI-SR is %4.4f\n',rho8');
fprintf('The PI-E is %4.4f\n',E8');
fprintf('The PI-EE is %4.4f\n',EE8');
fprintf('The ABC2-SR is %4.4f\n',rho9');
fprintf('The ABC2-E is %4.4f\n',E9');
fprintf('The ABC2-EE is %4.4f\n',EE9');
fprintf('The GA2-SR is %4.4f\n',rho10');
fprintf('The GA2-E is %4.4f\n',E10');
fprintf('The GA2-EE is %4.4f\n',EE10');