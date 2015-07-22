function DW=DebyeWaller(A, U_matrix,hkl)
a=A(:,1); b=A(:,2); c=A(:,3);
%Define the metric matrix according to p68 Giacovazzo
G=[a'*a a'*b a'*c; b'*a b'*b b'*c ; c'*a c'*b c'*c];
G_star=G^-1;   %using crystallographers definition to better follow Paltinus dissertation p 35
A_star=G_star*A; %p 144 Giacovazzo

M=2*pi^2*A_star*U_matrix*A_star;

DW=exp(-dot(hkl,M*hkl));

end