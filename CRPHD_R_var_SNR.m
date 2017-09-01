%%����֮���CRPHD��ͳ������,��������У����SNR�ı仯
clc;
close all;
clear all;

%% �źŲ���
N = 200;
n = 0:1:N-1;
% L = 90;

w = 0.3*pi;
V1 =1;
V2 = 2;
phi1 = 0*pi;
phi2 = 0*pi;
sn = V1*exp(1j*w*n - 1j*phi1) +  V2*exp(-1j*w*n -1j* phi2);

%% ��������
SNRV= -10 :5: 40;
M = 1000;
wAcu = w * ones(1,M);             %��ʵ���ź�
index = 1;
for SNR=SNRV
    for m = 1:M
        xn = awgn(sn,SNR);
        rnx = r(xn, N,N/2);
        rns = r(sn, N,N/2);
        RNX = length(rnx);
        RNS = length(rns);
        qn = rnx-rns;
        %% ��һ�ַ���
        %% ���Ƴ�����ֵ
        A1=real(AR(rnx,RNX));
        B1=real(BR(rnx,RNX));
        argumentc = (B1+sqrt(B1^2+8*A1^2))/(4*A1); 
        if (argumentc<-1)
            display('-1')
            argumentc = -1;
        end
        if (argumentc>1)
            display('1')
            argumentc = 1;
         end
         CRPHD_gu(m) = acos(argumentc);           
        %% ����ֵ
         A=real(AR(rns,RNS));
         B=real(BR(rns,RNS));
         kA = real(kexiA( rns,qn, RNS ));
         kB = real(kexiB( rns,qn, RNS ));
         k1 = 16*A*kA + 2*B*kB;
         x0 = 8*A^2+B^2;
         k2 = k1/(2*sqrt(x0));
         k3 = (-kA*sqrt(x0) + k2*A - B*kA + A*kB)/(4*A^2);
         y0 = (B*A + A*sqrt(x0))/(4*A^2);
         Fy(m) = acos(y0) - k3/sqrt(1-y0^2);
    end
    MSE_CRPHD(index) = 10*log10(mse(wAcu, CRPHD_gu));

    MSE_Fy(index) = 10*log10(mse(wAcu, Fy));
    index = index + 1;
    index
end
% %CRLB�ο��½�
% CRLB = 12 ./ ((4*pi*pi) * 10.^(SNRV/10) * N *(N^2-1));
% 
% plot(SNRV, 10*log10(CRLB), 'LineWidth', 2);

plot(SNRV,  MSE_CRPHD, '*', 'LineWidth', 2);
hold on;
plot(SNRV,  MSE_Fy, '-', 'LineWidth', 2);
%�����ʽ����
% axis([-10  40 -110 -20]);
legend('CRPHD','Theoretical');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
xlabel('\fontname{Times New Roman}SNR', 'FontWeight','bold');%����Times New Roman���Ӵ�
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;