% 第一种方法统计特性kexiB
function res = kexiB( s,q, N )
res = conj(s(N))*q(N) + conj(q(N))*s(N) - conj(s(N-1))*q(N-1) - conj(q(N-1))*s(N-1) - conj(s(2))*q(2) - conj(q(2))*s(2) + conj(s(1))*q(1) + conj(q(1))*s(1) ;
for n = 3:N
    res = res + 2*(conj(s(n))*q(n-2) + conj(q(n))*s(n-2));
end
end