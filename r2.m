% �ڶ��ַ����õ�������غ���
function res = r2( x, N, M,L)
 res = zeros(N-M+1,1);

    for k=0:N-M
        for n = 1:M-L
             for m = 0:L-1
                res(k+1) = res(k + 1) + x(n+m)*conj(x(n+m+k));
            end
        end
        res(k+1) = res(k + 1) / M;
    end
 