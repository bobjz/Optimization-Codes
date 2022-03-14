function [M] = makeM(a,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[~,d1]=size(a);
[~,d2]=size(b);
M=zeros(d1,d2);
for i=1:d1
    for j=1:d2
        M(i,j)=(a(i)-b(j))^2;
    end
end
end

