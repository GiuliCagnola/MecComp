close all; clear all;

% L1 y L2: extremos izquierdo y derecho respectivamente
% N: cantidad de intervalos (N+1 nodos)
% cb: condiciones de borde [tipo valor1 valor2]
  % tipo==1 -> Dirichet -- valor1 = Temp. impuesta
  % tipo==2 -> Neumann -- valor1 = q (flujo impuesto)
  % tipo==3 -> Robin -- valor1 = h (conveccion) -- valor2 = Temp. ext.
% model: struct con los parametros del modelo incluida la fuente
%        punto a punto
% et: esquema temporal [tipo val1 val2 val3]
  % tipo==0 -> Estacionario
  % tipo==1 -> Explicito
  % tipo==2 -> Implicito
  % Si tipo==1 o tipo==2: val1=maxIt, val2=tolError
  % Si tipo==2 val3=dt
%
% Para cada ejercicio se resuelve el estado estacionario, se compara con la
% solucion analitica (grafica + error por nodo) y se hace un estudio de
% convergencia: RMSE vs dx en escala log-log, cuya pendiente es el orden p
% del metodo (E ~ C*dx^p).

% ===================== CONFIGURACION =====================
ejercicios = 'abcdefg';          % ejercicios a correr, p.ej. 'c' o 'ad'
Nconv = [5 10 20 40 80 160];     % intervalos usados en el estudio de convergencia

% ===================== EJERCICIOS =====================
% Cada ejercicio devuelve: dominio, malla, CB, parametros, fuente G(x) y
% solucion analitica Tex(x) (ambas como funciones de x).
function ej = definirEjercicio(letra)
  ej.rho = 1; ej.cp = 1;
  switch letra
    case 'a'  # 2 T'' + 100 = 0, T(0)=10, T(1)=50
      ej.L1 = 0; ej.L2 = 1; ej.N = 100;
      ej.cb = [1 10 -1; 1 50 -1];
      ej.k = 2; ej.c = 0;
      ej.G   = @(x) 100*ones(size(x));
      ej.Tex = @(x) -25*x.^2 + 65*x + 10;

    case 'b'  # T'' - T = 0, T(0)=100, q(2)=0
      ej.L1 = 0; ej.L2 = 2; ej.N = 100;
      ej.cb = [1 100 -1; 2 0 -1];
      ej.k = 1; ej.c = 1;
      ej.G   = @(x) zeros(size(x));
      ej.Tex = @(x) (100*exp(-x).*(exp(2*x)+exp(4)))/(1+exp(4));

    case 'c'  # T'' + 100(x-3)^2 = 0, q(1)=2, T(5)=0
      ej.L1 = 1; ej.L2 = 5; ej.N = 100;
      ej.cb = [2 2 -1; 1 0 -1];
      ej.k = 1; ej.c = 0;
      ej.G   = @(x) 100*((x-3).^2);
      ej.Tex = @(x) (-25*x.^4 + 300*x.^3 - 1350*x.^2 + 1906*x + 2345)/3;

    case 'd'  # T'' - T + 50 = 0, T(0)=10, Robin en x=1 (h=0.2, Tinf=50)
      ej.L1 = 0; ej.L2 = 1; ej.N = 10;
      ej.cb = [1 10 -1; 3 0.2 50];
      ej.k = 1; ej.c = 1;
      ej.G   = @(x) 50*ones(size(x));
      # T = A e^-x + B e^x + 50, con A+B=-40 y la CB Robin -> B = 2A e^-2/3
      A = -40/(1 + 2*exp(-2)/3);  B = 2*A*exp(-2)/3;   # A=-36.6897, B=-3.3103
      ej.Tex = @(x) A*exp(-x) + B*exp(x) + 50;

    case 'e'  # 2 T'' + x^3 = 0, Robin en x=5 (h=2, Tinf=100), T(10)=50
      ej.L1 = 5; ej.L2 = 10; ej.N = 20;
      ej.cb = [3 2 100; 1 50 -1];
      ej.k = 2; ej.c = 0;
      ej.G   = @(x) x.^3;
      ej.Tex = @(x) -(x.^5)/40 + (1225*x)/3 - 4600/3;

    case 'f'  # 2 T'' - 2 T + 75 = 0, T(0)=0, Robin en x=1 (h=2, Tinf=10)
      ej.L1 = 0; ej.L2 = 1; ej.N = 100;
      ej.cb = [1 0 -1; 3 2 10];
      ej.rho = 1; ej.cp = 2;
      ej.k = 2; ej.c = 2;
      ej.G   = @(x) 75*ones(size(x));
      ej.Tex = @(x) (-5/4)*exp(-(x+1)).*(exp(x)-1).*(11*exp(x) + 11 - 30*exp(1));

    case 'g'  # 2 T'' + 2 T = 0, T(0)=50, q(1)=5
      ej.L1 = 0; ej.L2 = 1; ej.N = 100;
      ej.cb = [1 50 -1; 2 5 -1];
      ej.k = 2; ej.c = -2;
      ej.G   = @(x) zeros(size(x));
      # T = A sin(x) + 50 cos(x), con -k T'(1) = 5
      A = (50*sin(1) - 2.5)/cos(1);                    # A=73.2433
      ej.Tex = @(x) A*sin(x) + 50*cos(x);

    otherwise
      error('Ejercicio 2.%s no definido', letra);
  endswitch
endfunction

% Resuelve el estacionario con N intervalos y devuelve nodos, aproximada y exacta
function [xnode, T, Tex] = resolverEjercicio(ej, N)
  xnode = linspace(ej.L1, ej.L2, N+1)';
  model.rho = ej.rho; model.cp = ej.cp;
  model.k = ej.k; model.c = ej.c;
  model.G = ej.G(xnode);
  T = difFinitas(xnode, model, ej.cb, [0], false);
  Tex = ej.Tex(xnode);
endfunction

% ===================== RESOLUCION Y COMPARACION =====================
resumen = zeros(numel(ejercicios), 5);

for e = 1:numel(ejercicios)
  letra = ejercicios(e);
  ej = definirEjercicio(letra);

  % --- solucion con la malla del ejercicio ---
  [xnode, T, Tex] = resolverEjercicio(ej, ej.N);
  err  = T - Tex;
  rmse = sqrt(mean(err.^2));                 # raiz del error cuadratico medio
  emax = max(abs(err));                      # error maximo
  erel = norm(err)/norm(Tex);                # error relativo (norma 2)

  % --- estudio de convergencia ---
  dxC = (ej.L2 - ej.L1)./Nconv;
  rmseC = zeros(size(Nconv));
  for j = 1:numel(Nconv)
    [~, Tj, Texj] = resolverEjercicio(ej, Nconv(j));
    rmseC(j) = sqrt(mean((Tj - Texj).^2));
  endfor
  # si el error es de redondeo (solucion exacta polinomica de grado <= 3)
  # la pendiente no tiene sentido
  exactoMaq = all(rmseC < 1e-9*max(1, max(abs(Tex))));
  if exactoMaq
    orden = NaN;
  else
    p = polyfit(log(dxC), log(rmseC), 1);
    orden = p(1);
  endif

  resumen(e,:) = [ej.N, rmse, emax, erel, orden];

  % --- graficas ---
  figure('Name', sprintf('Ejercicio 2.%s', letra), 'Position', [100 100 1300 420]);

  subplot(1,3,1);
  xf = linspace(ej.L1, ej.L2, 400);
  plot(xf, ej.Tex(xf), 'r-', 'linewidth', 1.5); hold on;
  plot(xnode, T, 'bo', 'markersize', 4); hold off;
  grid on; xlabel('x'); ylabel('Temp');
  legend('Exacta', 'Dif. finitas');
  title(sprintf('Ejercicio 2.%s  (N = %d intervalos)', letra, ej.N));

  subplot(1,3,2);
  plot(xnode, abs(err), 'k.-');
  grid on; xlabel('x'); ylabel('|T_{ex} - T_{ap}|');
  title(sprintf('Error: RMSE = %.2e,  max = %.2e', rmse, emax));

  subplot(1,3,3);
  loglog(dxC, rmseC, 'bo-', 'linewidth', 1.5); hold on;
  if !exactoMaq
    # recta de referencia de pendiente 2 que pasa por el ultimo punto
    loglog(dxC, rmseC(end)*(dxC/dxC(end)).^2, 'r--');
    legend('RMSE', 'pendiente 2', 'location', 'northwest');
    title(sprintf('Convergencia: orden p = %.3f', orden));
  else
    title('Convergencia: exacta a precision de maquina');
  endif
  hold off; grid on; xlabel('\Delta x'); ylabel('RMSE');
endfor

% ===================== RESUMEN =====================
disp('---------------------------------------------------------------------');
printf('%-10s %6s %12s %12s %12s %8s\n', 'Ejercicio', 'N', 'RMSE', 'Error max', 'Error rel', 'Orden');
for e = 1:numel(ejercicios)
  printf('%-10s %6d %12.3e %12.3e %12.3e %8.3f\n', ['2.' ejercicios(e)], resumen(e,:));
endfor
disp('---------------------------------------------------------------------');
disp('Orden: pendiente de log(RMSE) vs log(dx). NaN = solucion exacta (error de redondeo).');

# Si se corre fuera de la GUI (p.ej. "octave archivo.m" desde VS Code),
# Octave cierra las figuras al terminar el script: esperamos a que el
# usuario las cierre.
if !isguirunning()
  disp('Cerrar las figuras para terminar...');
  for h = findobj('type', 'figure')'
    waitfor(h);
  endfor
endif
