%%������������е�����SNR�ı仯
clc;
close all;
clear all;

%% �źŲ���
N = 200;
n = 0:1:N-1;
% L = 95;    %�ڶ����ַ����ƴ�
l1=1;
l2=floor(N/3);%�����ַ����ֹ���
p=2;  %%�����������M��
q = 20;
w = 0.4*pi;
V1 = sqrt(2);
phi1 = 0*pi;
sn = V1*cos(w*n + phi1);
%% ��������
SNRV= -10 :2: 40;
M = 1000;
wAcu1 = w * ones(1,M);             %��ʵ���ź�
index = 1;
for SNR=SNRV
  
    for m = 1:M
      
        xn = awgn(sn,SNR);
%% ��һ�������
    rn1 = r1(xn, N,N/2);
    RN1 = length(rn1);
    A1=real(AC(rn1,RN1));
    B1=BC(rn1,RN1);
    argumentc1 =A1/(2*B1);
    if (argumentc1<-1)
        display('-1')
        argumentc1 = -1;
    end
    if (argumentc1>1)
        display('1')
        argumentc1 = 1;
     end
     CLS_gu1(m) = acos(argumentc1);
     %% �ڶ����ƴ������
%      rnX2 = r1(xn, N,N/2);
%         %% �ƴ����� 
%        for   kk = 1: N/2-L+1
%             rn2  = [ rnX2(kk:(kk+L-1))];
%             RN2 = length(rn2);
%             A2=real(AC(rn2,RN2));
%             B2=BC(rn2,RN2);
%             argumentc2 =A2/(2*B2);
%             if (argumentc2<-1)
%                 display('-1')
%                 argumentc2 = -1;
%             end
%             if (argumentc2>1)
%                 display('1')
%                 argumentc2 = 1;
%              end
%            CLS_gu20(kk) = acos(argumentc2);  
%        end
%        CLS_gu2(:,m) = CLS_gu20;
%    %% �������ƴ�����غ�����
%    rnX3 = r1(xn, N,N/2);
%    %% �ƴ�����
%    for   kk = 1: N/2-L+1
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
    %% �����������N-k
    for ii=1:N-1                  %%�������������
        Rn(ii)=r3(xn,N,ii);
    end
    rn4 = [Rn(l1:l2)];   %%��ȡ����������Լ���������Ӱ��
    RN4 = length(rn4);
    A4=real(AC( rn4, RN4));
    B4=real(BC(rn4, RN4));
    argumentc4 =A4/(2*B4);
    if (argumentc4<-1)
        display('-1')
        argumentc4 = -1;
    end
    if (argumentc4>1)
        display('1')
        argumentc4 = 1;
     end
     CLS_gu4(m) = acos(argumentc4);
     %% ��������ƫ�����
     for ii=1:N-1                  %%�������������
        Rn5(ii)=r4(xn,N,ii);
     end
        rn5 = [Rn5(p:q)];   %%��ȡ����������Լ���������Ӱ��
        RN5 = length(rn5);

        A5=real(AC( rn5, RN5));
        B5=real(BC(rn5, RN5));
        argumentc5 =A5/(2*B5);
        if (argumentc5<-1)
            display('-1')
            argumentc5 = -1;
        end
        if (argumentc5>1)
            display('1')
            argumentc5 = 1;
        end
     CLS_gu5(m) = acos(argumentc5);
       %% CRPHD��������
    
      %% ��һ�������
    rn6 = r1(xn, N,N/2);
    RN6 = length(rn6);
    A6=real(AR(rn6,RN6));
    B6=real(BR(rn6,RN6));
    argumentc6 = (B6+sqrt(B6^2+8*A6^2))/(4*A6); 
    if (argumentc6<-1)
        display('-1')
        argumentc6 = -1;
    end
    if (argumentc6>1)
        display('1')
        argumentc6 = 1;
     end
     CRPHD_gu6(m) = acos(argumentc6);
%      %% �������ƴ�����غ�����
%         rnX7 = r1(xn, N,N/2);
%         %% �ƴ��Ĳ���
%          for   kk = 1: N/2-L+1
%             rnx7 = [ rnX7(kk:(kk+L-1))];
%             RNX7 = length(rnx7);
%             for n = 0:RNX7-2
%                 Tn7(n+1) = rnx7(n+1)-rnx7(n+2);
%             end
%             TN7 = length(Tn7);
%             A7=real(AR(Tn7,TN7));
%             B7=real(BR(Tn7,TN7));
%             argumentc7 = (B7+sqrt(B7^2+8*A7^2))/(4*A7); 
%             if (argumentc7<-1)
%                 display('-1')
%                 argumentc7 = -1;
%             end
%             if (argumentc7>1)
%                 display('1')
%                 argumentc7 = 1;
%              end
%              CRPHD_gu70(kk) = acos(argumentc7);
%          end
%       CRPHD_gu7(:,m) = CRPHD_gu70;

    end
    MSE_CLS1(index) = 10*log10(mse(wAcu1, CLS_gu1));
%     MSE_CLS2(index) = 10*log10(mse(wAcu2, CLS_gu2));
%     MSE_CLS3(index) = 10*log10(mse(wAcu2, CLS_gu3));
    MSE_CLS4(index) = 10*log10(mse(wAcu1, CLS_gu4));
    MSE_CLS5(index) = 10*log10(mse(wAcu1, CLS_gu5));
    MSE_CLS6(index) = 10*log10(mse(wAcu1, CRPHD_gu6));
%     MSE_CLS7(index) = 10*log10(mse(wAcu2, CRPHD_gu7));
    index = index + 1;
    index
end
%CRLB�ο��½�
% CRLB = 12 ./ (4*(pi.^2)*10.^(SNRV/10) * N *(N^2-1));
CRLB = 12 ./ (10.^(SNRV/10) * N *(N^2-1));

plot(SNRV,  MSE_CLS1, '*', 'LineWidth', 2);
hold on;
% plot(SNRV,  MSE_CLS2, '-^', 'LineWidth', 2);
% plot(SNRV,  MSE_CLS3, '-^', 'LineWidth', 2);
plot(SNRV,  MSE_CLS4, '^', 'LineWidth', 2);
plot(SNRV,  MSE_CLS5, 'd', 'LineWidth', 2);
plot(SNRV,  MSE_CLS6, 'o', 'LineWidth', 2);
% plot(SNRV,  MSE_CLS7, '-^', 'LineWidth', 2);
plot(SNRV, 10*log10(CRLB),'k-', 'LineWidth', 2);
%�����ʽ����
set(gca, 'xtick', -10:5:40);
legend('CLS','Close-form expended','Multiple autocorrelation lags','CRPHD','CRLB');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
%legend('��һ�������','�ڶ���������ƴ�','������������ƴ�������','�����������N-k');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
xlabel('\fontname{Times New Roman}SNR', 'FontWeight','bold');%����Times New Roman���Ӵ�
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;
