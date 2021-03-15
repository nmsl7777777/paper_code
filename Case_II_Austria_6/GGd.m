function result=GGd(J,Z,M)
%J [M^2,1] => [M^2,Ni]
%Gd*J  [M^2,1] => [M^2,Ni]
Ni = size(J,2);
result = zeros(M^2,Ni);
for ii=1:Ni
temp1 = ifft2(fft2(Z).*fft2(reshape(J(:,ii),M,M),2*M-1,2*M-1));
temp2 = reshape(temp1(1:M,1:M),M^2,1);
result(:,ii) = temp2;
end