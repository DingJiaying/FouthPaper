%%����֮���CLS����������У���ͳ������
clc;
close all;
clear all;

%% �źŲ���
N = 200;

% L = 90;

wv = 0 : pi/25 : pi;   %��pi/25��pi-pi/25ÿ��pi/25ȡһ��
V1 =1;
V2 = 2;
phi1 = 0.0*pi;
phi2 = 0.0*pi;


%% ��������
SNR = 20;
M = 1000;
index = 1;
for w = wv 
    wAcu = w * ones(1,M);             %��ʵ���ź�
    for m = 1:M
        n = 0:1:N-1;
        sn = V1*exp(1j*w*n - 1j*phi1) +  V2*exp(-1j*w*n - 1j* phi2);
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
         CLS_gu1(m) = acos(argumentc);
        %% ���۵�ֵ
         C = real(AC(rns,RNS));
         D = 2*real(BC(rns,RNS)); 
         kC =  real(kexiC( rns,qn, RNS ));
         kD =  real(kexiD( rns,qn, RNS ));
         z0 = C/(D);
         k4 = (-C*kD+D*kC)/(D^2);
         Fz(m) = acos(z0) - k4/sqrt(1-z0^2);
%          %% �ƴ��Ĳ���
%        %% �������ƴ�����غ�����
%     rnX3 = r1(xn, N,N/2);
%     %% �ƴ�����
%     for   kk = 1: N/2-L+1
%         rn3  = [ rnX3(kk:(kk+L-1))];
%         RN3 = length(rn3);
%         for n = 0:RN3-2
%             Tn3(n+1) = rn3(n+1)-rn3(n+2);
%         end
%         TN3 = length(Tn3);
%         A3=real(AC(Tn3,TN3));
%         B3=BC(Tn3,TN3);
%         argumentc3 =A3/(2*B3);
%         if (argumentc3<-1)
%             display('-1')
%             argumentc3 = -1;
%         end
% 
%         if (argumentc3>1)
%             display('1')
%             argumentc3 = 1;
%         end
%          CLS_gu30(kk) = acos(argumentc3);
%     end
%     CLS_gu3(:,m) = CLS_gu30;
%     Lkk = length(CLS_gu30);
%     wAcu2 = w * ones(Lkk,M);             %��ʵ���ź�
    end
    MSE_CLS1(index) = 10*log10(mse(wAcu, CLS_gu1));
    MSE_Fz(index) = 10*log10(mse(wAcu, Fz));
%     MSE_CLS3(index) = 10*log10(mse(wAcu2, CLS_gu3));
    index = index + 1;
    index
end
plot(wv/pi,  MSE_CLS1, '*', 'LineWidth', 2);
hold on;
% plot(wv/pi,  MSE_CLS3, '-^', 'LineWidth', 2);
plot(wv/pi,  MSE_Fz, '-', 'LineWidth', 2);

%�����ʽ����
axis([0 1 -90 -60]);
legend('CLS','Theoretical');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
xlabel('\fontname{Times New Roman}\omega / \pi', 'FontWeight','bold');%����Times New Roman���Ӵ�
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;
