% 第一种方法统计特性kexiD
function res = kexiD( s,q, N )
res = 0;
for n = 3:N
    res = res + 2*(conj(s(n-1))*q(n-1) + conj(q(n-1))*s(n-1));
end
end