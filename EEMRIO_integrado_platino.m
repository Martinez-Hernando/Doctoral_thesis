As=readmatrix('As.txt');
n=length(As);
Ys=zeros(n,1);
Ys(111)=1;
Rs=readmatrix('Rs.txt');
I=eye(n-111);
Aio=readmatrix('Aio_zeros.txt');
int=I-Aio;
As(112:end,112:end)=int;
A_inv=inv(As);
Matriz=A_inv*Ys;
Matriz_rep= repmat(Matriz.', 18, 1);
Impacto_inicial=Rs*(A_inv*Ys);
impacto=Rs.*Matriz_rep;

%%RECICLADO
As_mod=As;
platino_total=abs(As(107,111)+As(108,111)+As(109,111)+As(110,111));
platino_rec=0.25*platino_total;
As_mod(109,111)=-platino_rec;
As_mod(107,111)=As(107,111)+(abs(As_mod(109,111))-abs(As(109,111)))*abs(As(107,111))/(abs(As(107,111)+As(108,111)));
As_mod(108,111)=As(108,111)+(abs(As_mod(109,111))-abs(As(109,111)))*abs(As(108,111))/(abs(As(107,111)+As(108,111)));
A_inv_mod=inv(As_mod);
Matriz5=A_inv_mod*Ys;
Matriz5_rep= repmat(Matriz5.', 18, 1);
Impacto4=Rs*(A_inv_mod*Ys);
impacto_rec_rep=Rs.*Matriz5_rep;
impacto_reciclado=Rs*(A_inv_mod*Ys);


%writematrix(Impacto_inicial,'Platino_metal.txt')
%writematrix(impacto,'Platino_metal_categorias.txt')
%writematrix(impacto_reciclado,'Platino_metal_rec.txt')
%writematrix(impacto_rec_rep,'Platino_metal_rec_categorias.txt')