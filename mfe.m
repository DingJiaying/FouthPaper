%  mean frequency error��ʽ ����Ϊmfe
function res = mfe(acu, arg, M)
res = 0;
for m = 1:M
   res = res + arg(m) - acu(m); 
end
res = res / M;
end