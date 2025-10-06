As=readmatrix('As_litio.txt');
%Yst=readmatrix('Ys.txt');
n=length(As);
Ys=zeros(n,1);
Ys(177)=1;
Rs=readmatrix('Rs_litio.txt');
I=eye(n-177);
Aio=readmatrix('Aio_litio_metal.txt');
int=I-Aio;
As(178:end,178:end)=int;
A_inv=inv(As);
Matriz=A_inv*Ys;
Matriz_rep= repmat(Matriz.', 18, 1);
impacto_inicial=Rs*(A_inv*Ys);
impacto=Rs.*Matriz_rep;

%%EN PORTUGAL
As_mod=As;
litio_total=abs(As(168,177)+As(169,177)+As(170,177)+As(171,177));
litio_portugal=0.1*litio_total;
As_mod(171,177)=-litio_portugal;
As_mod(168,177)=As(168,177)*0.9;
As_mod(169,177)=As(169,177)*0.9;
As_mod(170,177)=As(170,177)*0.9;

A_inv_mod=inv(As_mod);
Matriz5=A_inv_mod*Ys;
Matriz5_rep= repmat(Matriz5.', 18, 1);
impacto_port_rep=Rs.*Matriz5_rep;
impacto_portugal=Rs*(A_inv_mod*Ys);


%writematrix(impacto_inicial,'litio_metal.txt')
%writematrix(impacto,'litio_metal_categorias.txt')
%writematrix(impacto_portugal,'litio_metal_portugal.txt')
%writematrix(impacto_port_rep,'litio_metal_port_categorias.txt')