% 第一种方法kexiC
function res = kexiC( s,q, N )
res = 0;
for n = 3:N
   res = res + conj(s(n-1))*(q(n)+q(n-2))+conj(q(n-1))*(s(n)+s(n-2)); 
end
end