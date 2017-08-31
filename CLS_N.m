%%测试N对估计性能的影响
%%整理之后的CRPHD针对自相关序列
clc;
close all;
clear all;

%% 信号参数
N = 600;
n = 0:1:N-1;
L =100;     %%窗口长度
w = 0.3*pi;
V1 =1;
V2 = 2;
phi1 = 0.1*pi;
phi2 = 0.2*pi;
sn = V1*exp(1j*w*n - 1j* phi1) +  V2*exp(-1j*w*n - 1j*phi2);

SNR1 = -10;
SNR2 =0;
SNR3 =20;
SNR4 = 40;
SNR5 = 60;
%% 独立运行

M = 1000;
NN = 10:20:500;

index = 1;
for nn = NN
    wAcu = w * ones(1,M);             %真实的信号
    for m = 1:M
        xn1 = awgn(sn, SNR1);         %RPHD的含有加性高斯白噪声的信号
        xn2 = awgn(sn, SNR2); 
        xn3 = awgn(sn, SNR3); 
        xn4 = awgn(sn, SNR4); 
        xn5 = awgn(sn, SNR5); 
        
        rn1 = r(xn1, nn,nn/2);
        rn2 = r(xn2, nn,nn/2);
        rn3 = r(xn3, nn,nn/2);
        rn4 = r(xn4, nn,nn/2);
        rn5 = r(xn5, nn,nn/2);
        
        %CRPHD1的关键参数
        RN1 = length(rn1);
        A1=real(AC(rn1,RN1));
        B1=real(BC(rn1,RN1));
        argumentc1 = A1/(2*B1); 
        if (argumentc1<-1)
            display('-1')
            argumentc1 = -1;
        end
        if (argumentc1>1)
            display('1')
            argumentc1 = 1;
         end
         CRPHD_gu1(m) = acos(argumentc1);

         %CRPHD2的关键参数
        RN2 = length(rn2);
        A2=real(AC(rn2,RN2));
        B2=real(BC(rn2,RN2));
        argumentc2 = A2/(2*B2); 
        if (argumentc2<-1)
            display('-1')
            argumentc2 = -1;
        end
        if (argumentc2>1)
            display('1')
            argumentc2 = 1;
         end
         CRPHD_gu2(m) = acos(argumentc2);
         
         %CRPHD3的关键参数
        RN3 = length(rn3);
        A3=real(AC(rn3,RN3));
        B3=real(BC(rn3,RN3));
        argumentc3 = A3/(2*B3); 
        if (argumentc3<-1)
            display('-1')
            argumentc3 = -1;
        end
        if (argumentc3>1)
            display('1')
            argumentc3 = 1;
         end
         CRPHD_gu3(m) = acos(argumentc3);
         
         %CRPHD4的关键参数
        RN4 = length(rn4);
        A4=real(AC(rn4,RN4));
        B4=real(BC(rn4,RN4));
        argumentc4 = A4/(2*B4); 
        if (argumentc4<-1)
            display('-1')
            argumentc4 = -1;
        end
        if (argumentc4>1)
            display('1')
            argumentc4 = 1;
         end
         CRPHD_gu4(m) = acos(argumentc4);
         
         %CRPHD5的关键参数
        RN5 = length(rn5);
        A5=real(AC(rn5,RN5));
        B5=real(BC(rn5,RN5));
        argumentc5 = A5/(2*B5); 
        if (argumentc5<-1)
            display('-1')
            argumentc5 = -1;
        end
        if (argumentc5>1)
            display('1')
            argumentc5 = 1;
         end
         CRPHD_gu5(m) = acos(argumentc5);
         
    end
    MSE_CLS1(index) = 10*log10(mse(wAcu, CRPHD_gu1));
    MSE_CLS2(index) = 10*log10(mse(wAcu, CRPHD_gu2));
    MSE_CLS3(index) = 10*log10(mse(wAcu, CRPHD_gu3));
    MSE_CLS4(index) = 10*log10(mse(wAcu, CRPHD_gu4));
    MSE_CLS5(index) = 10*log10(mse(wAcu, CRPHD_gu5));
    index
    index = index + 1;
end
plot(NN, MSE_CLS1, '-^', 'LineWidth', 2);
hold on;
plot(NN, MSE_CLS2, '-o', 'LineWidth', 2);
plot(NN, MSE_CLS3, '-*', 'LineWidth', 2);
plot(NN, MSE_CLS4, '-d', 'LineWidth', 2);
plot(NN, MSE_CLS5, '->', 'LineWidth', 2);
%仿真格式部分
legend('SNR= -10dB','SNR=0dB','SNR=20dB','SNR=40dB','SNR=60dB');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例
xlabel('\fontname{Times New Roman}N', 'FontWeight','bold');%字体Times New Roman，加粗
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;