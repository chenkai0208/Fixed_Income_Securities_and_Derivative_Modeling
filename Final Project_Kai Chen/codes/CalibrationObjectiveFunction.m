% This is just the set up for the optimization problem. To solve the
% system, we recommend readers use Matlab function @fminsearch dedicated
% for nonlinear optimization by simply putting in the command window:
% [Lambda, f] = fminsearch(@CalibrationObjectiveFunction, Lambda0);

function f = CalibrationObjectiveFunction(Lambda)
Sig = [0.227  0.23   0.221 0.209  0.196  0.186  0.176 0.169  0.163  0.159;  % test data
       0.224  0.215  0.205 0.194  0.183  0.174  0.167 0.162  0.158  0.154;
       0.209  0.201  0.19  0.18   0.17   0.163  0.158 0.155  0.152  0.15;
       0.195  0.187  0.177 0.168  0.16   0.155  0.151 0.148  0.147  0.145;
       0.182  0.174  0.165 0.158  0.151  0.148  0.145 0.143  0.142  0.14;
       0.1746 0.1674 0.159 0.1524 0.1462 0.1436 0.141 0.1394 0.1384 0.1368;
       0.1672 0.1608 0.153 0.1468 0.1414 0.1392 0.137 0.1358 0.1348 0.1336;
       0.1598 0.1542 0.147 0.1412 0.1366 0.1348 0.133 0.1322 0.1312 0.1304;
       0.1524 0.1476 0.141 0.1356 0.1318 0.1304 0.129 0.1286 0.1276 0.1272;
       0.145  0.141  0.135 0.13   0.127  0.1260 0.125 0.125  0.124  0.124];

%Lambda = [1 2 3 4 5 6 7 8 9 10];  % input as an initial guess

Today='25-Jan-2005';
B=[0.9774658 0.9509789 0.9219838 0.8911017 0.8591725 0.8264399 0.7930540 0.7597502 0.7262834 0.6944457 0.6645450 0.6349818 0.6068399 0.5792752 0.5523236 0.5273147 0.5030900 0.4795796 0.4567881 0.4346590];
VectorOfDates = ['25-Jan-2006'; '25-Jan-2007'; '25-Jan-2008'; '26-Jan-2009'; '25-Jan-2010'; '25-Jan-2011'; '25-Jan-2012'; '25-Jan-2013'; '27-Jan-2014'; '26-Jan-2015'; '25-Jan-2016'; '25-Jan-2017'; '25-Jan-2018'; '25-Jan-2019'; '27-Jan-2020'; '25-Jan-2021'; '25-Jan-2022'; '25-Jan-2023'; '25-Jan-2024'; '27-Jan-2025'];
T_Num = datenum(VectorOfDates);
m = 10; % Number of swaption maturities(ex.10)
M = 20; % Number of swaption maturities plus number of swaption underlyings 
R = []; % Setting zeros for matrix [R] as initial values
for i=1:m
    for j=i+1:M-m+i
        for k=i+1:j
            R(i,j,k)=(B(k-1)-B(k))/(B(i)-B(j));
        end
    end
end

VCV = [];
for k=1:m
    VCV(k,k)= yearfrac(Today,T_Num(k))*Sig(k,1)^2/Lambda(k); 
end
s=1;
for i=1:m
    for j=i+1:m 
        Sum=0;
        for l=i+j-2*s+1:j+1 
            for k=i+j-2*s+1:j+1
                SumTemp=R(i+j-2*s,j+1,k)*R(i+j-2*s,j+1,1)*VCV(k-1,l-1);
                Sum = Sum + SumTemp; 
            end
        end
        VCV(i+j-2*s,j) = (yearfrac(Today,T_Num(i+j-2*s))*Sig(i+j-2*s,i+1)^2 - Lambda(i+j-2*s)*(Sum-2*R(i+j-2*s,j+1,i+j-2*s+1)*VCV(i+j-2*s,j)*R(i+j-2*s,j+1,j+1)))/(2*Lambda(i+j-2*s)*R(i+j-2*s,j+1,i+j-2*s+1)*R(i+j-2*s,j+1,j+1));
        VCV(j,i+j-2*s)=VCV(i+j-2*s,j); 
    end
    s=s+1;
end

[E,X]=eig(VCV);
L = diag(X);

for i=1:m
    if L(i)<0
        L_check(i)=1;
    else
        L_check(i)=0;
    end
end

for i=1:m
    if L_check(i)==0
        for j=1:m
            E_sqrL(j,i)=E(j,i)*sqrt(L(i));
        end
    else
        for j=1:m
            E_sqrL(j,i)=0;
        end
    end
end
VCV_M=E_sqrL*E_sqrL';

Sig_theo=[];
for k=1:m
    for N=k+1:m+1
        Sum=0;
        for l=k+1:N
            for i=k+1:N
                SumTemp=R(k,N,i)*VCV_M(i-1,l-1)*R(k,N,l);
                Sum = Sum+SumTemp;
            end
        end
        Sig_theo(k,N-k)=sqrt(Sum*Lambda(k)/yearfrac(Today,T_Num(k)));
    end
end

save dat Sig_theo

RSME = 0; 
for i = 1:m
    for j=1:m-i+1
        RSME_Temp=(Sig_theo(i, j)-Sig(i, j))^2; 
        RSME = RSME+RSME_Temp;
    end
end
f = RSME; % function f will be used as a minimization function









