% =========================================================
% TEST: verificacion contra los items a, b, c de la guia
% =========================================================


fprintf('\n=== ITEM a: dominio [0,1], k=2, c=0, G=100, T(0)=10, T(1)=50 ===\n')
xnode = linspace(0, 1, 3);
model.k=2; model.c=0; model.rho=0; model.cp=1; model.G=100;
cb = [1, 10, -1;   % Dirichlet T=10 en x=0
      1, 50, -1];  % Dirichlet T=50 en x=1
T = elemFinitos(xnode, model, cb, 0);
Tanal = @(x) -25*x.^2 + 65*x + 10;
fprintf('  x      T_MEF   T_anal  error\n')
for i=1:length(xnode)
    fprintf('  %.2f   %6.3f  %6.3f  %6.4f\n', xnode(i), T(i), Tanal(xnode(i)), abs(T(i)-Tanal(xnode(i))))
end

fprintf('\n=== ITEM b: dominio [0,2], k=1, c=1, G=0, T(0)=100, q(2)=0 ===\n')
xnode = linspace(0, 2, 3);
model.k=1; model.c=1; model.rho=0; model.cp=1; model.G=0;
cb = [1, 100, -1;  % Dirichlet T=100
      2, 0,   -1]; % Neumann q=0
T = elemFinitos(xnode, model, cb, 0);
Tanal_b = @(x) 100*cosh(2-x)/cosh(2);   % -T''+T=0, T(0)=100, T'(2)=0
fprintf('  x      T_MEF   T_anal  error\n')
for i=1:length(xnode)
    fprintf('  %.2f   %6.3f  %6.3f  %6.4f\n', xnode(i), T(i), Tanal_b(xnode(i)), abs(T(i)-Tanal_b(xnode(i))))
end

fprintf('\n=== ITEM d: dominio [0,1], k=1, c=1, rho*cp=1, G=1, T(0)=10, Robin en x=1 ===\n')
xnode = linspace(0, 1, 3);
model.k=1; model.c=1; model.rho=1; model.cp=1; model.G=1;
cb = [1, 10,  -1;   % Dirichlet T=10
      3, 0.2, 50];  % Robin hcv=0.2, Tinf=50
T_FE = elemFinitos(xnode, model, cb, 1);
T_BE = elemFinitos(xnode, model, cb, 2);
T_CN = elemFinitos(xnode, model, cb, 3);
% estacionario: -T''+T=1 -> T = 1 + A*e^x + B*e^-x
%   T(0)=10  y  -T'(1) = hcv*(T(1)-Tinf)
hcv = 0.2; Tinf = 50;
AB = [1, 1; (1+hcv)*exp(1), (hcv-1)*exp(-1)] \ [10-1; hcv*(Tinf-1)];
Tanal_d = @(x) 1 + AB(1)*exp(x) + AB(2)*exp(-x);
fprintf('  x      FE      BE      CN      T_anal\n')
for i=1:length(xnode)
    fprintf('  %.2f   %6.3f  %6.3f  %6.3f  %6.3f\n', xnode(i), T_FE(i), T_BE(i), T_CN(i), Tanal_d(xnode(i)))
end

fprintf('\nTodos los tests completados.\n')