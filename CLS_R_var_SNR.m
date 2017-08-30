%%����֮���CLS����������У���ͳ������
clc;
close all;
clear all;

%% �źŲ���
N = 200;
n = 0:1:N-1;


w = 0.3*pi;
V1 = 1;
V2 = 2;
phi1 = 0.1*pi;
phi2 = 0.2*pi;
sn = V1*exp(1j*w*n - 1j*phi1) +  V2*exp(-1j*w*n - 1j*phi2);

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
        %% ���Ƶ�ֵ
        C1 = real(AC(rnx, RNX));
        D1 = BC(rnx, RNX);
        argumentc =C1/(2*D1);
        if (argumentc<-1)
            display('-1')
            argumentc = -1;
        end
        if (argumentc>1)
            display('1')
            argumentc = 1;
         end
         CLS_gu(m) = acos(argumentc);

          %% ���۵�ֵ
         C = real(AC(rns,RNS));
         D = 2*real(BC(rns,RNS)); 
         kC =  real(kexiC( rns,qn, RNS ));
         kD =  real(kexiD( rns,qn, RNS ));
         z0 = C/(D);
         k4 = (-C*kD+D*kC)/(D^2);
         Fz(m) = acos(z0) - k4/sqrt(1-z0^2);
    end
    MSE_CLS(index) = 10*log10(mse(wAcu, CLS_gu));
     MSE_Fz(index) = 10*log10(mse(wAcu, Fz));

    %% ѭ������
         index = index + 1;
    index
end
plot(SNRV,  MSE_CLS, '*', 'LineWidth', 2);
hold on;
% plot(SNRV,  MSE_CLS1, '-^', 'LineWidth', 2);
plot(SNRV,  MSE_Fz, '-', 'LineWidth', 2);
%�����ʽ����
legend('CLS','Theoretical');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
xlabel('\fontname{Times New Roman}SNR', 'FontWeight','bold');%����Times New Roman���Ӵ�
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;
