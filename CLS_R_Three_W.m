%%三种方法随着W的变化，针对RCLS
clc;
close all;
clear all;

%% 信号参数
N = 200;
% L = 90;    %第二三种方法移窗
l1=1;
l2=floor(N/3);%第四种方法粗估计
p=2;  %%第五种自相关M阶
q = 20;
wv = 0 : pi/25 : pi;   %从pi/25到pi-pi/25每隔pi/25取一个
V1 = sqrt(2);
phi1 = 0*pi;

%% 独立运行
SNR = 20;
M = 1000;

index = 1;
for w = wv 
    wAcu1 = w * ones(1,M);             %真实的信号
    for m = 1:M
        n = 0:1:N-1;
        sn = V1*cos(w*n + phi1);
        xn = awgn(sn,SNR);
%% 第一种自相关
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
  
%     %% 第三种移窗自相关后项差分
%     rnX3 = r1(xn, N,N/2);
%     %% 移窗部分
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
%     wAcu2 = w * ones(Lkk,M);             %真实的信号
    %% 第四种自相关N-k
    for ii=1:N-1                  %%产生自相关序列
        Rn(ii)=r3(xn,N,ii);
    end
    rn4 = [Rn(l1:l2)];   %%截取自相关序列以减弱噪声的影响
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
     %% 第五种无偏自相关
     for ii=1:N-1                  %%产生自相关序列
        Rn5(ii)=r4(xn,N,ii);
     end
        rn5 = [Rn5(p:q)];   %%截取自相关序列以减弱噪声的影响
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
     %% CRPHD的两条线
    
      %% 第一种自相关
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
%      %% 第三种移窗自相关后项差分
%         rnX7 = r1(xn, N,N/2);
%         %% 移窗的部分
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
%     MSE_CLS3(index) = 10*log10(mse(wAcu2, CLS_gu3));
     MSE_CLS4(index) = 10*log10(mse(wAcu1, CLS_gu4));
     MSE_CLS5(index) = 10*log10(mse(wAcu1, CLS_gu5));
      MSE_CLS6(index) = 10*log10(mse(wAcu1, CRPHD_gu6));
%       MSE_CLS7(index) = 10*log10(mse(wAcu2, CRPHD_gu7));
    index = index + 1;
    index
end
%CRLB参考下界
CRLBtemp = 0;
w1 = 0*pi:0.01*pi:0.99*pi;
for k = 0:N-1
%     CRLBtemp = CRLBtemp + (k^2*(sin(w1*k+phi1).^2));
%     CRLBtemp = CRLBtemp + (k*2*pi*sin(w1*k+phi1).^2);
    CRLBtemp = CRLBtemp + (k*k/8*pi*sin(w1*k+phi1).^2);
end

CRLB = 1 ./ (10^(SNR/10)*CRLBtemp);


plot(wv/pi,  MSE_CLS1, '*', 'LineWidth', 2);
hold on;
% plot(wv/pi,  MSE_CLS3, '-^', 'LineWidth', 2);
plot(wv/pi,  MSE_CLS4, '^', 'LineWidth', 2);
plot(wv/pi,  MSE_CLS5, 'd', 'LineWidth', 2);
plot(wv/pi,  MSE_CLS6, 'o', 'LineWidth', 2);
% plot(wv/pi,  MSE_CLS7, '-^', 'LineWidth', 2);
plot(w1/pi, 10*log10(CRLB),'k-' ,'LineWidth', 2);
%仿真格式部分
% axis([0 1 -95 -65]);
%legend('第一种自相关','第二种自相关移窗','第三种自相关移窗后项差分','第四种自相关N-k');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例

legend('CLS','Close-form expended','Multiple autocorrelation lags','CRPHD','CRLB');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例
xlabel('\fontname{Times New Roman}\omega / \pi', 'FontWeight','bold');%字体Times New Roman，加粗
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;
