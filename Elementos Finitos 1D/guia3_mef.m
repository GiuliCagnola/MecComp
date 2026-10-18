close all; clear all; more off;
% =========================================================================
% GUIA DE TRABAJOS PRACTICOS N 3 - Metodo de Elementos Finitos 1D
%
%   rho*cp dT/dt = k T'' - c T + G     en [L1, L2]
%
% Para cada item:
%  1) Estacionario con elementos lineales y cuadraticos, con la misma
%     cantidad de nodos, comparado con la solucion analitica.
%  2) Convergencia: error en norma L2 contra la analitica al refinar.
%     Teoria: lineales ~ h^2, cuadraticos ~ h^3.
%  3) Items con rho*cp > 0 (e, f, g): transitorio desde T(x,0) = 0 con
%     Forward Euler, Backward Euler y Crank-Nicolson.
%
% El Ejercicio 2 de la guia es la funcion elemFinitos.m (en esta misma
% carpeta).
% =========================================================================

items  = 'abcdefg';          % items a resolver
nNodos = 9;                  % nodos para la comparacion (8 lineales / 4 cuadraticos)
neConv = [2 4 8 16 32 64];   % elementos para el estudio de convergencia
neTrans = 8;                 % elementos lineales para el transitorio

% ----------------------------------------------------------------------
% Datos de la tabla. cb = [tipo v1 v2]: [1 T -1] Dirichlet,
% [2 q -1] Neumann (q saliente), [3 h Tinf] Robin
% ----------------------------------------------------------------------
function it = definirItem(letra)
  switch letra
    case 'a'
      it.L = [0 1];  it.rhocp = 0; it.k = 2; it.c = 0; it.G = 100;
      it.cb = [1 10 -1; 1 50 -1];
      it.Tex = @(x) -25*x.^2 + 65*x + 10;
    case 'b'
      it.L = [0 2];  it.rhocp = 0; it.k = 1; it.c = 1; it.G = 0;
      it.cb = [1 100 -1; 2 0 -1];
      it.Tex = @(x) 100*exp(-x).*(exp(2*x) + exp(4))/(1 + exp(4));
    case 'c'
      it.L = [1 5];  it.rhocp = 0; it.k = 1; it.c = 0; it.G = @(x) 100*(x-3).^2;
      it.cb = [2 2 -1; 1 0 -1];
      it.Tex = @(x) (-25*x.^4 + 300*x.^3 - 1350*x.^2 + 1906*x + 2345)/3;
    case 'd'
      it.L = [0 1];  it.rhocp = 0; it.k = 1; it.c = 1; it.G = 50;
      it.cb = [1 10 -1; 3 0.2 50];
      A = -40/(1 + 2*exp(-2)/3);  B = 2*A*exp(-2)/3;     % -36.6897, -3.3103
      it.Tex = @(x) A*exp(-x) + B*exp(x) + 50;
    case 'e'
      it.L = [5 10]; it.rhocp = 1; it.k = 2; it.c = 0; it.G = @(x) x.^3;
      it.cb = [3 2 100; 1 50 -1];
      it.Tex = @(x) -x.^5/40 + 1225*x/3 - 4600/3;
    case 'f'
      it.L = [0 1];  it.rhocp = 2; it.k = 2; it.c = 2; it.G = 75;
      it.cb = [1 0 -1; 3 2 10];
      it.Tex = @(x) -5/4*exp(-(x+1)).*(exp(x) - 1).*(11*exp(x) + 11 - 30*exp(1));
    case 'g'
      it.L = [0 1];  it.rhocp = 1; it.k = 2; it.c = -2; it.G = 0;
      it.cb = [1 50 -1; 2 5 -1];
      A = (50*sin(1) - 2.5)/cos(1);                       % 73.2433
      it.Tex = @(x) A*sin(x) + 50*cos(x);
  endswitch
endfunction

function model = armarModelo(it, orden)
  model.k = it.k;  model.c = it.c;  model.G = it.G;
  model.rho = it.rhocp;  model.cp = 1;  model.orden = orden;
endfunction

% Solucion MEF evaluada en puntos xf (interpolando con las funciones de forma)
function Tf = interpolarMEF(xnode, T, orden, xf)
  xnode = xnode(:);  T = T(:);  Tf = zeros(size(xf));
  for i = 1 : numel(xf)
    e = find(xnode(1:orden:end-orden) <= xf(i) + 1e-12, 1, 'last');
    nod = (e-1)*orden + (1:orden+1);
    xe = xnode(nod);
    xi = 2*(xf(i) - xe(1))/(xe(end) - xe(1)) - 1;   % nodo medio centrado
    if orden == 1
      N = [(1-xi)/2, (1+xi)/2];
    else
      N = [xi*(xi-1)/2, 1-xi^2, xi*(xi+1)/2];
    end
    Tf(i) = N*T(nod);
  end
endfunction

% Error en norma L2: sqrt( integral (T_h - T_ex)^2 dx ), Gauss de 5 puntos
function err = errorL2(xnode, T, orden, Tex)
  xg = [-0.9061798459 -0.5384693101 0 0.5384693101 0.9061798459];
  wg = [0.2369268851 0.4786286705 0.5688888889 0.4786286705 0.2369268851];
  xv = xnode(1:orden:end);  s = 0;
  for e = 1 : numel(xv)-1
    xq = (xv(e)+xv(e+1))/2 + (xv(e+1)-xv(e))/2*xg;
    d = interpolarMEF(xnode, T, orden, xq) - Tex(xq);
    s = s + (xv(e+1)-xv(e))/2 * sum(wg.*d.^2);
  end
  err = sqrt(s);
endfunction

% ======================================================================
for letra = items
  it = definirItem(letra);
  L1 = it.L(1);  L2 = it.L(2);
  xf = linspace(L1, L2, 400);
  printf('\n================ Item %s ================\n', letra);

  esTrans = it.rhocp > 0;
  figure('Name', sprintf('Ejercicio 1.%s', letra), 'Position', [50 50 1450 420]);
  nsub = 2 + esTrans;

  % ---------------- 1) Estacionario: lineales vs cuadraticos ----------------
  xnode = linspace(L1, L2, nNodos);
  Tlin  = elemFinitos(xnode, armarModelo(it, 1), it.cb, 0);
  Tcuad = elemFinitos(xnode, armarModelo(it, 2), it.cb, 0);
  Tan   = it.Tex(xnode(:));
  printf('Estacionario, %d nodos (%d lineales / %d cuadraticos):\n', nNodos, nNodos-1, (nNodos-1)/2);
  printf('%9s %12s %12s %12s %11s %11s\n', 'x', 'lineal', 'cuadratico', 'analitica', 'err lin', 'err cuad');
  printf('%9.4f %12.4f %12.4f %12.4f %11.2e %11.2e\n', [xnode(:) Tlin Tcuad Tan abs(Tlin-Tan) abs(Tcuad-Tan)]');

  subplot(1, nsub, 1);
  plot(xf, it.Tex(xf), 'k-', 'linewidth', 2); hold on;
  plot(xf, interpolarMEF(xnode, Tlin, 1, xf), 'r--', 'linewidth', 1.2);
  plot(xf, interpolarMEF(xnode, Tcuad, 2, xf), 'b-.', 'linewidth', 1.2);
  plot(xnode, Tlin, 'ro', xnode, Tcuad, 'b+', 'markersize', 7);
  hold off; grid on; xlabel('x'); ylabel('T');
  legend('analitica', sprintf('%d lineales', nNodos-1), sprintf('%d cuadraticos', (nNodos-1)/2), ...
         'location', 'southoutside', 'orientation', 'horizontal');
  title(sprintf('Item %s: estacionario', letra));

  % ---------------- 2) Convergencia en norma L2 ----------------
  subplot(1, nsub, 2);
  he = (L2 - L1)./neConv;
  estilo = {'ro-', 'bs-'};  leyenda = {};  p = [NaN NaN];
  for orden = 1 : 2
    e = zeros(size(neConv));
    for j = 1 : numel(neConv)
      xn = linspace(L1, L2, orden*neConv(j) + 1);
      Tj = elemFinitos(xn, armarModelo(it, orden), it.cb, 0);
      e(j) = errorL2(xn, Tj, orden, it.Tex);
    end
    printf('Error L2 %-11s %s', {'lineal:', 'cuadratico:'}{orden}, num2str(e, '%10.2e'));
    ok = e > 1e-10*max(abs(Tan));          % se descarta el error de redondeo
    if nnz(ok) < 2
      % la solucion exacta esta en el espacio de aproximacion (item a con
      % cuadraticos): el error es solo de redondeo
      printf('   -> exacta (error de redondeo)\n');
      continue;
    end
    q = polyfit(log(he(ok)), log(e(ok)), 1);  p(orden) = q(1);
    printf('   -> p = %.2f\n', p(orden));
    loglog(he(ok), e(ok), estilo{orden}, 'linewidth', 1.5); hold on;
    leyenda{end+1} = sprintf('%s (p = %.2f)', {'lineales', 'cuadraticos'}{orden}, p(orden));
  end
  hold off; grid on; xlabel('h_e'); ylabel('||T - T_h||_{L2}');
  legend(leyenda, 'location', 'northwest');
  if isnan(p(2))
    title('Convergencia (cuadraticos: exacta)');
  else
    title('Convergencia');
  end

  % ---------------- 3) Transitorio: FE, BE y CN ----------------
  if esTrans
    xnode = linspace(L1, L2, neTrans + 1);
    model = armarModelo(it, 1);
    [~, iM] = min(abs(xnode - (L1+L2)/2));          % punto de control: x medio
    Test = elemFinitos(xnode, model, it.cb, 0);     % estacionario MEF
    % dt critico de Forward Euler (se lo reporta elemFinitos con el esquema 1)
    [~, ~, tFE] = elemFinitos(xnode, model, it.cb, 1);
    dtFE = tFE(2);
    printf('Transitorio (%d lineales), T(x,0) = 0, punto x = %g:\n', neTrans, xnode(iM));
    printf('  %-16s %9s %7s %9s %14s\n', 'esquema', 'dt', 'pasos', 't final', 'max|T-T_est|');

    subplot(1, nsub, 3); hold on;
    casos = {1, dtFE, 'r-', 'FE'; 2, dtFE, 'b--', 'BE'; 3, dtFE, 'g-.', 'CN'; ...
             2, 20*dtFE, 'bo', 'BE'; 3, 20*dtFE, 'gs', 'CN'};
    leyenda = {};
    for m = 1 : rows(casos)
      [T, H, t] = elemFinitos(xnode, model, it.cb, [casos{m,1} casos{m,2}]);
      plot(t, H(iM,:), casos{m,3}, 'linewidth', 1.3, 'markersize', 5);
      leyenda{end+1} = sprintf('%s, dt = %.3g', casos{m,4}, casos{m,2});
      printf('  %-16s %9.4g %7d %9.3f %14.2e\n', casos{m,4}, casos{m,2}, numel(t)-1, t(end), max(abs(T-Test)));
    end
    xl = [0 1.5*tFE(end)];               % hasta un poco despues de que FE llega
    xlim(xl);
    plot(xl, it.Tex(xnode(iM))*[1 1], 'k:', 'linewidth', 1.5);
    leyenda{end+1} = 'estacionario analitico';
    hold off; grid on; xlabel('t'); ylabel(sprintf('T(%g, t)', xnode(iM)));
    legend(leyenda, 'location', 'southeast');
    title(sprintf('Transitorio (dt critico FE = %.3g)', dtFE/0.9));
  end
  drawnow;
end

% Si se corre fuera de la GUI (p.ej. "octave archivo.m" desde VS Code),
% Octave cierra las figuras al terminar el script: esperamos a que el
% usuario las cierre.
if !isguirunning()
  disp('Cerrar las figuras para terminar...');
  for h = findobj('type', 'figure')'
    waitfor(h);
  endfor
endif
