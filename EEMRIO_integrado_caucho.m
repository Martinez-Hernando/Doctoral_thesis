As=readmatrix('As_caucho.txt');
n=length(As);
Ys=zeros(n,1);
Ys(78)=1;
Rs=readmatrix('Rs_caucho.txt');
I=eye(n-78);
Aio=readmatrix('Aio_caucho.txt');
int=I-Aio;
As(79:end,79:end)=int;
As_1=As;
As_1(79:end,29:51)=0;
A_inv=inv(As_1);
Matriz=A_inv*Ys;
Matriz_rep=repmat(Matriz.', 18, 1);
impacto_inicial=Rs*(A_inv*Ys);
impacto=Rs.*Matriz_rep;

%%EN ESPAÑA
Aio_esp=readmatrix('Aio_caucho_esp.txt');
int_mod=I-Aio_esp;
As_mod=As;
As_mod(79:end,79:end)=int_mod;
caucho_total=abs(As(74,78)+As(75,78)+As(76,78)+As(77,78));
caucho_spain=0.1*caucho_total;
As_mod(77,78)=-caucho_spain;
As_mod(74,78)=As(74,78)*0.9;
As_mod(75,78)=As(75,78)*0.9;
As_mod(76,78)=As(76,78)*0.9;

A_inv_mod=inv(As_mod);
Matriz5=A_inv_mod*Ys;
Matriz5_rep= repmat(Matriz5.', 18, 1);
impacto_esp_rep=Rs.*Matriz5_rep;
impacto_esp=Rs*(A_inv_mod*Ys);


%writematrix(impacto_inicial,'litio_metal.txt')
%writematrix(impacto,'litio_metal_categorias.txt')
%writematrix(impacto_portugal,'litio_metal_portugal.txt')
%writematrix(impacto_port_rep,'litio_metal_port_categorias.txt')