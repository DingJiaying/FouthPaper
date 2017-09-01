%%����֮���CRPHD��ͳ������,��������У�MFE
clc;
close all;
clear all;

%% �źŲ���
N = 200;

L = 90;

wv = pi/25 : pi/25 : pi-pi/25;   %��pi/25��pi-pi/25ÿ��pi/25ȡһ��
V1 =1;
V2 = 2;
phi1 = 0.1*pi;
phi2 = 0.2*pi;


%% ��������
SNR = 10;
M = 100;

index = 1;
for w = wv 
    wAcu = w * ones(1,M);             %��ʵ���ź�
    for m = 1:M
        n = 0:1:N-1;
        sn = V1*exp(1j*w*n + phi1) +  V2*exp(-1j*w*n + phi2);
        xn = awgn(sn,SNR);
        rnx = r(xn, N,N/2);
        rns = r(sn, N,N/2);
        RNX = length(rnx);
        RNS = length(rns);
        qn = rnx-rns;
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
         %% �������ƴ�����غ�����
        rnX3 = r1(xn, N,N/2);
        %% �ƴ��Ĳ���
         for   kk = 1: N/2-L+1
            rnx3 = [ rnX3(kk:(kk+L-1))];
            RNX3 = length(rnx3);
            for n = 0:RNX3-2
                Tn3(n+1) = rnx3(n+1)-rnx3(n+2);
            end
            TN3 = length(Tn3);
            A3=real(AR(Tn3,TN3));
            B3=real(BR(Tn3,TN3));
            argumentc3 = (B3+sqrt(B3^2+8*A3^2))/(4*A3); 
            if (argumentc3<-1)
                display('-1')
                argumentc3 = -1;
            end
            if (argumentc3>1)
                display('1')
                argumentc3 = 1;
             end
             CRPHD_gu20(kk) = acos(argumentc3);
         end
      CRPHD_gu2(:,m) = CRPHD_gu20;
       Lkk = length(CRPHD_gu20);
    end
    MSE_CRPHD(index) = mfe(wAcu, CRPHD_gu, M);
    MSE_CRPHD2(index) = mfe(wAcu, CRPHD_gu2, M);

    index = index + 1;
    index
end
plot(wv/pi,  MSE_CRPHD, '^', 'LineWidth', 2);
hold on;
plot(wv/pi,  MSE_CRPHD2, 'o', 'LineWidth', 2);

line([0,1], [0,0],'LineWidth', 2,'color','k');
%�����ʽ����
set(gca, 'xtick', 0:0.1:1);
% set(gca, 'ytick', -0.00004:0.000005:0.00004);
legend('Eq.10','Eq.13');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
xlabel('\fontname{Times New Roman}\omega / \pi', 'FontWeight','bold');%����Times New Roman���Ӵ�
ylabel('\fontname{Times New Roman}MFE (dB)', 'FontWeight','bold');
grid on;
hold off;