close all; clear all; more off;
% =========================================================================
% COMPARATIVA 1D: DIFERENCIAS FINITAS vs ELEMENTOS FINITOS
% Ejercicio 1 de las guias N 1 (DF) y N 3 (MEF): mismos items a-g.
%
%   rho*cp dT/dt = k T'' - c T + G     en [L1, L2]
%
% Para cada item se resuelve con los MISMOS nodos usando:
%   DF       : difFinitas  (stencil de 3 puntos, CB con nodo ficticio)
%   MEF lin  : elemFinitos con elementos lineales
%   MEF cuad : elemFinitos con elementos cuadraticos
% y se compara:
%   1) la solucion estacionaria contra la analitica,
%   2) el error en cada nodo,
%   3) la convergencia del error nodal (RMSE en los nodos) al refinar,
%   4) items con rho*cp > 0 (e, f, g): el transitorio con Euler implicito y
%      el mismo dt. DF usa rho*cp*I (masa "concentrada") y MEF la matriz de
%      capacidad consistente C, que acopla cada nodo con sus vecinos.
% =========================================================================
raiz = fileparts(mfilename('fullpath'));          % carpeta Codigos
addpath(fullfile(raiz, 'Diferencias Finitas 1D'));   % difFinitas, dirchlet
addpath(fullfile(raiz, 'Elementos Finitos 1D'));     % elemFinitos

items  = 'abcdefg';             % items a resolver ('p' = problema propio, ver abajo)
nNodos = 9;                     % nodos para la comparacion directa
nConv  = 2.^(2:7) + 1;          % nodos para el estudio de convergencia (impares)

% ----------------------------------------------------------------------
% Datos de la tabla de las guias. cb = [tipo v1 v2]:
%   [1 T -1] Dirichlet, [2 q -1] Neumann (q saliente), [3 h Tinf] Robin
% ----------------------------------------------------------------------
function it = definirItem(letra)
  switch letra
    case 'a'
      it.L = [0 1];  it.rhocp = 0; it.k = 2; it.c = 0;  it.G = @(x) 100*ones(size(x));
      it.cb = [1 10 -1; 1 50 -1];
      it.Tex = @(x) -25*x.^2 + 65*x + 10;
    case 'b'
      it.L = [0 2];  it.rhocp = 0; it.k = 1; it.c = 1;  it.G = @(x) zeros(size(x));
      it.cb = [1 100 -1; 2 0 -1];
      it.Tex = @(x) 100*exp(-x).*(exp(2*x) + exp(4))/(1 + exp(4));
    case 'c'
      it.L = [1 5];  it.rhocp = 0; it.k = 1; it.c = 0;  it.G = @(x) 100*(x-3).^2;
      it.cb = [2 2 -1; 1 0 -1];
      it.Tex = @(x) (-25*x.^4 + 300*x.^3 - 1350*x.^2 + 1906*x + 2345)/3;
    case 'd'
      it.L = [0 1];  it.rhocp = 0; it.k = 1; it.c = 1;  it.G = @(x) 50*ones(size(x));
      it.cb = [1 10 -1; 3 0.2 50];
      A = -40/(1 + 2*exp(-2)/3);  B = 2*A*exp(-2)/3;     % -36.6897, -3.3103
      it.Tex = @(x) A*exp(-x) + B*exp(x) + 50;
    case 'e'
      it.L = [5 10]; it.rhocp = 1; it.k = 2; it.c = 0;  it.G = @(x) x.^3;
      it.cb = [3 2 100; 1 50 -1];
      it.Tex = @(x) -x.^5/40 + 1225*x/3 - 4600/3;
      it.dt = 0.1;   it.tf = 20;                       % transitorio
    case 'f'
      it.L = [0 1];  it.rhocp = 2; it.k = 2; it.c = 2;  it.G = @(x) 75*ones(size(x));
      it.cb = [1 0 -1; 3 2 10];
      it.Tex = @(x) -5/4*exp(-(x+1)).*(exp(x) - 1).*(11*exp(x) + 11 - 30*exp(1));
      it.dt = 0.01;  it.tf = 2;
    case 'g'
      it.L = [0 1];  it.rhocp = 1; it.k = 2; it.c = -2; it.G = @(x) zeros(size(x));
      it.cb = [1 50 -1; 2 5 -1];
      A = (50*sin(1) - 2.5)/cos(1);                       % 73.2433
      it.Tex = @(x) A*sin(x) + 50*cos(x);
      it.dt = 0.01;  it.tf = 4;

    case 'p'
      % ================= PROBLEMA PROPIO (p.ej. el del parcial) =================
      % Completar los datos y poner items = 'p' arriba.
      it.L     = [];      % [L1 L2]               ej: [0 2]
      it.rhocp = [];      % rho*cp (0 = solo estacionario)
      it.k     = [];      % conductividad (constante)
      it.c     = [];      % reaccion (constante; negativa = fuente proporcional a T)
      it.G     = [];      % fuente: @(x) ... vectorizada (.* ./ .^); ej: @(x) 20*sin(x) + x.^2
                          %   si es constante alcanza con @(x) 50
      it.cb    = [];      % 2x3, fila 1 = izquierda, fila 2 = derecha:
                          %   [1 T -1] Dirichlet | [2 q -1] Neumann (q SALIENTE) | [3 h Tinf] Robin
                          %   ej: [1 15 -1; 3 4 30]
      it.Tex   = [];      % solucion analitica @(x) ... (OPCIONAL: si queda vacia se
                          %   compara contra una referencia con malla muy fina)
      it.dt    = [];      % solo si rho*cp > 0: paso de tiempo     ej: 0.01
      it.tf    = [];      % solo si rho*cp > 0: tiempo final       ej: 3
    otherwise
      error('Item "%s" no definido.', letra);
  endswitch
endfunction

% Valida los datos del item y completa lo que falta:
%  - G y Tex devuelven un valor por punto aunque se escriban constantes
%  - sin solucion analitica: referencia con MEF cuadratico de 1025 nodos
function it = completarItem(it, letra)
  req = {'L', 'rhocp', 'k', 'c', 'G', 'cb'};
  if isfield(it, 'rhocp') && !isempty(it.rhocp) && it.rhocp > 0
    req = [req, {'dt', 'tf'}];
  end
  for r = req
    if !isfield(it, r{1}) || isempty(it.(r{1}))
      error('Item %s: falta definir it.%s en definirItem.', letra, r{1});
    end
  end
  if numel(it.L) ~= 2 || it.L(2) <= it.L(1)
    error('Item %s: it.L debe ser [L1 L2] con L1 < L2.', letra);
  end
  if !isequal(size(it.cb), [2 3])
    error('Item %s: it.cb debe ser una matriz de 2x3.', letra);
  end
  if !isa(it.G, 'function_handle')
    it.G = @(x) it.G;               % G dada como numero
  end
  G0 = it.G;
  it.G = @(x) G0(x) + zeros(size(x));          % un valor por punto

  if isfield(it, 'Tex') && !isempty(it.Tex)
    T0 = it.Tex;
    it.Tex = @(x) T0(x) + zeros(size(x));
    it.nombreRef = 'analitica';
  else
    % referencia: MEF cuadratico con 1025 nodos. Las mallas de la comparacion
    % (2^m + 1 nodos) son subconjuntos de esta, asi que en sus nodos la
    % interpolacion devuelve exactamente el valor de referencia.
    xr = linspace(it.L(1), it.L(2), 1025)';
    mef.rho = 0;  mef.cp = 1;  mef.k = it.k;  mef.c = it.c;  mef.G = it.G;  mef.orden = 2;
    Tr = elemFinitos(xr, mef, it.cb, 0);
    it.Tex = @(x) reshape(interp1(xr, Tr, x(:), 'spline'), size(x));
    it.nombreRef = 'referencia (MEF cuad., 1025 nodos)';
    printf('Item %s: sin solucion analitica -> se compara contra una referencia con 1025 nodos.\n', letra);
  end
endfunction

% Resuelve el estacionario con los tres metodos sobre los nodos x
function [Tdf, Tlin, Tcuad] = resolverTres(it, x)
  x = x(:);
  mdf.rho = 1;  mdf.cp = 1;  mdf.k = it.k;  mdf.c = it.c;  mdf.G = it.G(x);
  Tdf = difFinitas(x, mdf, it.cb, [0], false);
  mef.rho = it.rhocp;  mef.cp = 1;  mef.k = it.k;  mef.c = it.c;  mef.G = it.G;
  mef.orden = 1;  Tlin  = elemFinitos(x, mef, it.cb, 0);
  mef.orden = 2;  Tcuad = elemFinitos(x, mef, it.cb, 0);
endfunction

nombres = {'DF', 'MEF lineal', 'MEF cuadratico'};
estilos = {'rs', 'bo', 'g^'};
lineas  = {'rs-', 'bo-', 'g^-'};
tabla = zeros(numel(items), 6);   % [errores nodales maximos (N=9) , ordenes]

for ii = 1 : numel(items)
  letra = items(ii);
  it = completarItem(definirItem(letra), letra);
  L1 = it.L(1);  L2 = it.L(2);
  esTrans = it.rhocp > 0;
  if esTrans
    figure('Name', sprintf('Item %s', letra), 'Position', [50 50 1200 820]);
    sp = @(i) subplot(2, 2, i);
  else
    figure('Name', sprintf('Item %s', letra), 'Position', [50 50 1450 430]);
    sp = @(i) subplot(1, 3, i);
  end

  % ---------------- 1) Solucion con los mismos nodos ----------------
  x = linspace(L1, L2, nNodos)';
  [Tdf, Tlin, Tcuad] = resolverTres(it, x);
  Tan = it.Tex(x);
  sols = {Tdf, Tlin, Tcuad};

  sp(1);
  xf = linspace(L1, L2, 400);
  plot(xf, it.Tex(xf), 'k-', 'linewidth', 2); hold on;
  for m = 1 : 3
    plot(x, sols{m}, estilos{m}, 'markersize', 8 - m, 'linewidth', 1.2);
  end
  hold off; grid on; xlabel('x'); ylabel('T');
  legend([it.nombreRef, nombres], 'location', 'southoutside', 'orientation', 'horizontal');
  title(sprintf('Item %s: estacionario con %d nodos', letra, nNodos));

  % ---------------- 2) Error en cada nodo ----------------
  sp(2);
  piso = 1e-14*max(abs(Tan));      % para poder graficar errores nulos en log
  for m = 1 : 3
    err = abs(sols{m} - Tan);
    semilogy(x, max(err, piso), lineas{m}, 'linewidth', 1.2, 'markersize', 6); hold on;
    tabla(ii, m) = max(err);
  end
  hold off; grid on; xlabel('x'); ylabel('|T - T_{ex}|');
  legend(nombres, 'location', 'southoutside', 'orientation', 'horizontal');
  title('Error nodal');

  % ---------------- 3) Convergencia del error nodal ----------------
  sp(3);
  h = (L2 - L1)./(nConv - 1);
  rmse = zeros(3, numel(nConv));
  for j = 1 : numel(nConv)
    xj = linspace(L1, L2, nConv(j))';
    [a, b, c] = resolverTres(it, xj);
    Tj = it.Tex(xj);
    rmse(:,j) = sqrt(mean(([a b c] - Tj).^2))';
  end
  leyenda = {};
  for m = 1 : 3
    ok = rmse(m,:) > 1e-10*max(abs(Tan));        % se descarta el error de redondeo
    if nnz(ok) >= 2
      q = polyfit(log(h(ok)), log(rmse(m,ok)), 1);  p = q(1);
      loglog(h(ok), rmse(m,ok), lineas{m}, 'linewidth', 1.5); hold on;
      leyenda{end+1} = sprintf('%s (p = %.2f)', nombres{m}, p);
    else
      p = NaN;                                    % exacta en los nodos
      leyenda{end+1} = sprintf('%s: exacta en nodos', nombres{m});
      loglog(NaN, NaN, lineas{m}); hold on;       % solo para la leyenda
    end
    tabla(ii, 3+m) = p;
  end
  hold off; grid on; xlabel('h'); ylabel('RMSE en los nodos');
  legend(leyenda, 'location', 'southeast');
  title('Convergencia (error nodal)');

  % ---------------- 4) Transitorio: Euler implicito, mismo dt ----------------
  if esTrans
    nP = round(it.tf/it.dt);
    % DF: tol = 0 -> hace exactamente nP pasos
    mdf.rho = it.rhocp;  mdf.cp = 1;  mdf.k = it.k;  mdf.c = it.c;  mdf.G = it.G(x);
    Hdf = difFinitas(x, mdf, it.cb, [2 nP 0 it.dt], false);
    % MEF lineal y cuadratico: Backward Euler hasta tf
    mef.rho = it.rhocp;  mef.cp = 1;  mef.k = it.k;  mef.c = it.c;  mef.G = it.G;
    mef.orden = 1;  [~, Hlin, t]  = elemFinitos(x, mef, it.cb, [2 it.dt it.tf]);
    mef.orden = 2;  [~, Hcuad]    = elemFinitos(x, mef, it.cb, [2 it.dt it.tf]);
    iM = (nNodos + 1)/2;                          % nodo central
    sp(4);
    plot(t, Hdf(iM,:), 'r-', t, Hlin(iM,:), 'b--', t, Hcuad(iM,:), 'g-.', 'linewidth', 1.5);
    hold on;
    plot([0 it.tf], it.Tex(x(iM))*[1 1], 'k:', 'linewidth', 1.5);
    hold off; grid on; xlabel('t'); ylabel(sprintf('T(%g, t)', x(iM)));
    legend([nombres, ['estacionario: ' it.nombreRef]], 'location', 'southeast');
    title(sprintf('Transitorio: Euler implicito, dt = %g, T(x,0) = 0', it.dt));
    dif = max(abs(Hdf(iM,:) - Hlin(iM,:)));
    printf('Item %s, transitorio en x = %g: max |T_DF - T_MEFlin| a lo largo del tiempo = %.3f\n', ...
           letra, x(iM), dif);
  end
  drawnow;
end

% ======================= Resumen =======================
printf('\n%s\n', repmat('-', 1, 92));
printf('%-5s | %-40s | %-40s\n', 'Item', sprintf('Error nodal maximo (%d nodos)', nNodos), 'Orden de convergencia (error nodal)');
printf('%-5s | %12s %12s %12s  | %12s %12s %12s\n', '', 'DF', 'MEF lin', 'MEF cuad', 'DF', 'MEF lin', 'MEF cuad');
printf('%s\n', repmat('-', 1, 92));
for ii = 1 : numel(items)
  printf('%-5s | %12.2e %12.2e %12.2e  | %12.2f %12.2f %12.2f\n', items(ii), tabla(ii,:));
end
printf('%s\n', repmat('-', 1, 92));
printf('NaN = exacta en los nodos (el error es solo de redondeo).\n');

% Si se corre fuera de la GUI (p.ej. "octave archivo.m" desde VS Code),
% Octave cierra las figuras al terminar el script: esperamos a que el
% usuario las cierre.
if !isguirunning()
  disp('Cerrar las figuras para terminar...');
  for h = findobj('type', 'figure')'
    waitfor(h);
  endfor
endif
