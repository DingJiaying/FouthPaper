%%����֮���CLS
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
    A=real(AC(xn,N));
    B=BC(xn,N);
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
 
% %%����һ�ַ������þ������ʽ
% SNR = 40;
% M = 100;
% wAcu = w0 * ones(1,M);             %��ʵ���ź�
% for m = 1:M
%     xn = awgn(sn,SNR);
%     xn = xn.';
%     for kk  = 1:N-L-1
%         v_2 = [xn(kk+2:(kk+L+1))];  %s(n+1)
%         v_0 = [xn(kk+1:(kk+L))];  %s(n)
%         v_1 = [xn(kk:(kk+L-1))];  %s(n-1)
%         A = [v_0'*(v_2+ v_1)];
%         B = [(v_0'*v_0)];
%         argumentc =pinv(2*B)*real(A);
%     end
%     if (argumentc<-1)
%         display('-1')
%         argumentc = -1;
%     end
%     if (argumentc>1)
%         display('1')
%         argumentc = 1;
%      end
%      CLS_gu(m) = acos(argumentc);
% end
% MSE_CLS = 10*log10(mse(wAcu, CLS_gu));