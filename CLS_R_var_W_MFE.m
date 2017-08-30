%%整理之后的CLS，自相关序列，的统计特性MFE
clc;
close all;
clear all;

%% 信号参数
N = 200;
% L = 95;

wv = pi/25 : pi/25 : pi-pi/25;   %从pi/25到pi-pi/25每隔pi/25取一个
V1 =1;
V2 = 2;
phi1 = 0*pi;
phi2 = 0*pi;


%% 独立运行
SNR = 10;
M = 10;

index = 1;
for w = wv 
    wAcu = w * ones(1,M);             %真实的信号
    for m = 1:M
        n = 0:1:N-1;
        sn = V1*exp(1j*w*n - 1j*phi1) +  V2*exp(-1j*w*n - 1j* phi2);
        xn = awgn(sn,SNR);
        rnx = r(xn, N,N/2);
        rns = r(sn, N,N/2);
        RNX = length(rnx);
        RNS = length(rns);
        qn = rnx-rns;
        %% 估计的值
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
%          %% 移窗的部分
%        %% 第三种移窗自相关后项差分
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
%          CLS_gu20(kk) = acos(argumentc3);

 %% CRPHD估计出来的值
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

    end
%     CLS_gu2(:,m) = CLS_gu20;
%     Lkk = length(CLS_gu20);
%     wAcu2 = w * ones(Lkk,M);             %真实的信号
%     end
    MSE_CLS(index) = mfe(wAcu, CLS_gu, M);
    MSE_CRPHD(index) = mfe(wAcu, CRPHD_gu, M);
    index = index + 1;
    index
end
plot(wv/pi,  MSE_CLS, 'd', 'LineWidth', 2);
hold on;
plot(wv/pi,  MSE_CRPHD,'o', 'LineWidth', 2);
line([0,1], [0,0],'LineWidth', 2,'color','k');
%仿真格式部分
set(gca, 'xtick', 0:0.1:1);
set(gca, 'ytick', -0.02:0.0002:0.02);
legend('CLS','CRPHD');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例
xlabel('\fontname{Times New Roman}\omega / \pi', 'FontWeight','bold');%字体Times New Roman，加粗
ylabel('\fontname{Times New Roman}MFE (dB)', 'FontWeight','bold');
grid on;
hold off;
