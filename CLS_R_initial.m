%%����֮��Ķ���������е�CLS
clc;
close all;
clear all;

%% �źŲ���
N = 200;
n = 0:1:N-1;
L =100;     %%���ڳ���
w0 = 0.3*pi;
V1 =1;
V2 = 2;
phi1 = 0.1*pi;
phi2 = 0.2*pi;
sn = V1*exp(1j*w0*n + phi1) +  V2*exp(-1j*w0*n + phi2);

%% ��������
SNR = 20;
M = 100;
wAcu = w0 * ones(1,M);             %��ʵ���ź�
for m = 1:M
    xn = awgn(sn,SNR);
    
    rn = r(xn, N,N/2);
    RN = length(rn);
    A=real(AC(rn,RN));
    B=BC(rn,RN);
    argumentc =A/(2*B);
    if (argumentc<-1)
        display('-1')
        argumentc = -1;
    end
    if (argumentc>1)
        display('1')
        argumentc = 1;
     end
     CLS_gu(m) = acos(argumentc);

end
MSE_CLS = 10*log10(mse(wAcu, CLS_gu))