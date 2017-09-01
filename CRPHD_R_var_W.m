%%整理之后的CRPHD的统计特性,自相关序列，MFE
clc;
close all;
clear all;

%% 信号参数
N = 200;

% L = 90;

wv =0 : pi/25 : pi;   %从pi/25到pi-pi/25每隔pi/25取一个
V1 =1;
V2 = 2;
phi1 = 0*pi;
phi2 = 0*pi;


%% 独立运行
SNR = 20;
M = 1000;

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
        %% 估计出来的值
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
%          %% 第三种移窗自相关后项差分
%         rnX3 = r1(xn, N,N/2);
%         %% 移窗的部分
%          for   kk = 1: N/2-L+1
%             rnx3 = [ rnX3(kk:(kk+L-1))];
%             RNX3 = length(rnx3);
%             for n = 0:RNX3-2
%                 Tn3(n+1) = rnx3(n+1)-rnx3(n+2);
%             end
%             TN3 = length(Tn3);
%             A3=real(AR(Tn3,TN3));
%             B3=real(BR(Tn3,TN3));
%             argumentc3 = (B3+sqrt(B3^2+8*A3^2))/(4*A3); 
%             if (argumentc3<-1)
%                 display('-1')
%                 argumentc3 = -1;
%             end
%             if (argumentc3>1)
%                 display('1')
%                 argumentc3 = 1;
%              end
%              CRPHD_gu20(kk) = acos(argumentc3);
%          end
%       CRPHD_gu2(:,m) = CRPHD_gu20;
%        Lkk = length(CRPHD_gu20);
%        wAcu2 = w * ones(Lkk,M);             %真实的信号
        %% 理论值
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
%     MSE_CRPHD2(index) = 10*log10(mse(wAcu2, CRPHD_gu2));
    MSE_Fy(index) = 10*log10(mse(wAcu, Fy));

    index = index + 1;
    index
end
plot(wv/pi,  MSE_CRPHD, '*', 'LineWidth', 2);
hold on;
% plot(wv/pi,  MSE_CRPHD2, '-^', 'LineWidth', 2);
plot(wv/pi,  MSE_Fy, '-', 'LineWidth', 2);

%仿真格式部分
axis([0 1 -90 -60]);
legend('CRPHD','Theoretical');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例
xlabel('\fontname{Times New Roman}\omega / \pi', 'FontWeight','bold');%字体Times New Roman，加粗
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;