%%整理之后的对自相关序列的CLS
clc;
close all;
clear all;

%% 信号参数
N = 200;
n = 0:1:N-1;
L =100;     %%窗口长度
w0 = 0.3*pi;
V1 =1;
V2 = 2;
phi1 = 0.1*pi;
phi2 = 0.2*pi;
sn = V1*exp(1j*w0*n + phi1) +  V2*exp(-1j*w0*n + phi2);

%% 独立运行
SNR = 20;
M = 100;
wAcu = w0 * ones(1,M);             %真实的信号
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