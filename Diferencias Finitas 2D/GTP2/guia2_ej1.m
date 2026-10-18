close all; clear all; more off;
% =========================================================================
% GUIA DE TRABAJOS PRACTICOS N 2 - Diferencias Finitas 2D - Ejercicio 1
%
%   rho*cp dT/dt = k (T_xx + T_yy) - c T + G   en [x1,x2] x [y1,y2]
%
% Para cada item:
%  1) Estado estacionario en mallas cada vez mas finas (n = 2, 4, 8, 16
%     intervalos por unidad de longitud) para:
%       - malla uniforme, CB de 2do orden (nodo ficticio)
%       - malla uniforme, CB de 1er orden (diferencia de 2 puntos)
%       - malla no uniforme (mas fina cerca de los bordes), CB de 2do orden
%     No hay solucion analitica: el orden se estima comparando mallas
%     sucesivas. Si e_h ~ C h^p, la diferencia entre la solucion con h y con
%     h/2 (evaluada en los nodos de la malla mas gruesa) tambien es ~ h^p, y
%       p = log2( dif(h, h/2) / dif(h/2, h/4) )
%  2) Items con rho*cp > 0 (e y f): transitorio con Forward Euler y Backward
%     Euler desde T(x,y,0) = 0, comparado con el estacionario.
% =========================================================================
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

items  = 'abcdef';       % items a resolver
nRef   = [2 4 8 16];     % intervalos por unidad de longitud
nGraf  = 8;              % malla usada para graficar y para el transitorio
alfaNU = 0.5;            % intensidad de la no uniformidad (0 = uniforme)

% ----------------------------------------------------------------------
% Datos de cada item. cb = [tipo v1 v2]:
%   [1 T 0] Dirichlet, [2 q 0] Neumann (q saliente), [3 h phi_inf] Robin
% ----------------------------------------------------------------------
function it = definirItem(letra)
  switch letra
    case 'a'
      it.ext = [0 0 1 1];  it.rhocp = 0; it.k = 2; it.c = 0;
      it.G = @(x,y) 300*ones(size(x));
      it.S = [1 20 0];  it.N = [1 20 0];  it.W = [1 20 0];  it.E = [1 20 0];
    case 'b'
      it.ext = [0 0 1 2];  it.rhocp = 0; it.k = 1; it.c = 4;
      it.G = @(x,y) 300*ones(size(x));
      it.S = [1 20 0];  it.N = [1 0 0];   it.W = [1 20 0];  it.E = [1 0 0];
    case 'c'
      it.ext = [0 0 4 3];  it.rhocp = 0; it.k = 1; it.c = 0;
      it.G = @(x,y) 2*x.^2 + 3*y.^3;
      it.S = [2 0 0];   it.N = [1 50 0];  it.W = [1 50 0];  it.E = [2 0 0];
    case 'd'
      it.ext = [1 1 3 3];  it.rhocp = 0; it.k = 1; it.c = 4;
      it.G = @(x,y) zeros(size(x));
      it.S = [1 25 0];  it.N = [3 2 50];  it.W = [2 0 0];   it.E = [2 10 0];
    case 'e'
      it.ext = [0 0 3 5];  it.rhocp = 1; it.k = 2; it.c = 0;
      it.G = @(x,y) 10*x.^3;
      it.S = [3 2 100]; it.N = [1 50 0];  it.W = [1 0 0];   it.E = [1 25 0];
      it.dtBE = [0.01 0.1];
    case 'f'
      it.ext = [0 0 2 2];  it.rhocp = 2; it.k = 2; it.c = 2;
      it.G = @(x,y) 100*ones(size(x));
      it.S = [1 0 0];   it.N = [3 2 10];  it.W = [3 2 10];  it.E = [1 100 0];
      it.dtBE = [0.01 0.1];
  endswitch
endfunction

% Arma malla, condiciones de borde y modelo de un item
function [xnode, icone, DIR, NEU, ROB, model, nx, ny] = armarProblema(it, n, alfa)
  x1 = it.ext(1); y1 = it.ext(2); x2 = it.ext(3); y2 = it.ext(4);
  nx = round(n*(x2-x1));  ny = round(n*(y2-y1));
  [xnode, icone] = malla_rect(x1, y1, x2, y2, nx, ny, alfa);
  s(1) = struct('dir', 1, 'pos', y1, 'rango', [x1 x2], 'cb', it.S);
  s(2) = struct('dir', 3, 'pos', y2, 'rango', [x1 x2], 'cb', it.N);
  s(3) = struct('dir', 4, 'pos', x1, 'rango', [y1 y2], 'cb', it.W);
  s(4) = struct('dir', 2, 'pos', x2, 'rango', [y1 y2], 'cb', it.E);
  [DIR, NEU, ROB] = cb_segmentos(xnode, icone, s);
  nn = rows(xnode);
  model.nnodes = nn;
  model.k = it.k*ones(nn,1);
  model.c = it.c*ones(nn,1);
  model.G = it.G(xnode(:,1), xnode(:,2));
  model.rho = it.rhocp;  model.cp = 1;
  model.ts = 2;  model.PHI_n = zeros(nn,1);
endfunction

% ----------------------------------------------------------------------
configs = struct('nombre', {'uniforme, CB 2do orden', 'uniforme, CB 1er orden', ...
                            'no uniforme, CB 2do orden'}, ...
                 'alfa', {0, 0, alfaNU}, 'orden', {2, 1, 2}, ...
                 'estilo', {'bo-', 'rs-', 'g^-'});

for letra = items
  it = definirItem(letra);
  tieneNR = any([it.S(1) it.N(1) it.W(1) it.E(1)] > 1);   % hay Neumann o Robin?
  xc = mean(it.ext([1 3]));  yc = mean(it.ext([2 4]));    % punto de control

  printf('\n================ Item %s ================\n', letra);
  printf('Dominio [%g,%g]x[%g,%g], rho*cp=%g, k=%g, c=%g\n', it.ext([1 3 2 4]), it.rhocp, it.k, it.c);

  esTrans = it.rhocp > 0;
  figure('Name', sprintf('Ejercicio 1.%s', letra), 'Position', [50 50 1400 430]);
  nsub = 2 + esTrans;

  % ---------------- 1) Estacionario y orden de convergencia ----------------
  subplot(1, nsub, 2); hold on;
  leyenda = {};
  printf('%-28s', 'T en el centro, n =');  printf('%10d', nRef);  printf('   orden p\n');
  for cfg = configs
    if cfg.orden == 1 && !tieneNR
      continue;     % solo Dirichlet: el orden de la CB no interviene
    end
    Tc = zeros(size(nRef));  sub = {};
    for j = 1 : numel(nRef)
      [xnode, icone, DIR, NEU, ROB, model, nx, ny] = armarProblema(it, nRef(j), cfg.alfa);
      model.cb_orden = cfg.orden;
      T = resolver(xnode, icone, DIR, NEU, ROB, model);
      [~, iC] = min((xnode(:,1)-xc).^2 + (xnode(:,2)-yc).^2);
      Tc(j) = T(iC);
      r = nRef(j)/nRef(1);                       % nodos de la malla mas gruesa
      Tg = reshape(T, nx+1, ny+1);
      sub{j} = Tg(1:r:end, 1:r:end);
      if cfg.alfa == 0 && cfg.orden == 2 && nRef(j) == nGraf
        Tgraf = T; xgraf = xnode; igraf = icone;
      end
    end
    dif = zeros(1, numel(nRef)-1);
    for j = 1 : numel(nRef)-1
      dif(j) = sqrt(mean((sub{j}(:) - sub{j+1}(:)).^2));
    end
    p = log2(dif(1:end-1)./dif(2:end));
    printf('%-28s', cfg.nombre);  printf('%10.4f', Tc);  printf('   %s\n', num2str(p, '%6.2f'));
    loglog(1./nRef(1:end-1), dif, cfg.estilo, 'linewidth', 1.5);
    leyenda{end+1} = sprintf('%s (p = %.2f)', cfg.nombre, p(end));
  end
  hh = 1./nRef(1:end-1);
  loglog(hh, dif(end)*(hh/hh(end)).^2, 'k--');
  leyenda{end+1} = 'pendiente 2';
  hold off; grid on; set(gca, 'xscale', 'log', 'yscale', 'log');
  xlabel('h = 1/n'); ylabel('RMS(T_h - T_{h/2})');
  legend(leyenda, 'location', 'northwest');
  title('Convergencia (estacionario)');

  subplot(1, nsub, 1);
  graficar_campo(xgraf, igraf, Tgraf, sprintf('Item %s: T estacionaria (n = %d)', letra, nGraf));

  % ---------------- 2) Transitorio: Forward Euler y Backward Euler ----------------
  if esTrans
    [xnode, icone, DIR, NEU, ROB, model] = armarProblema(it, nGraf, 0);
    [~, iC] = min((xnode(:,1)-xc).^2 + (xnode(:,2)-yc).^2);
    model.maxit = 200000;  model.tol = 1e-8;
    Test = Tgraf;

    subplot(1, nsub, 3); hold on;
    model.ts = 0;                                     % Forward Euler
    [Tfe, ~, dtFE] = resolver(xnode, icone, DIR, NEU, ROB, model);
    plot((0:columns(Tfe)-1)*dtFE, Tfe(iC,:), 'b-', 'linewidth', 2);
    leyenda = {sprintf('FE, dt = %.4g', dtFE)};
    printf('Transitorio (n = %d), T(x,y,0) = 0:\n', nGraf);
    printf('  Forward Euler  dt = %-8.4g pasos = %-6d t final = %-7.3f max|T - T_est| = %.1e\n', ...
           dtFE, columns(Tfe)-1, (columns(Tfe)-1)*dtFE, max(abs(Tfe(:,end)-Test)));
    estilos = {'r--', 'g-.'};
    for m = 1 : numel(it.dtBE)                        % Backward Euler
      model.ts = 1;  model.dt = it.dtBE(m);
      Tbe = resolver(xnode, icone, DIR, NEU, ROB, model);
      plot((0:columns(Tbe)-1)*model.dt, Tbe(iC,:), estilos{m}, 'linewidth', 1.5);
      leyenda{end+1} = sprintf('BE, dt = %g', model.dt);
      printf('  Backward Euler dt = %-8.4g pasos = %-6d t final = %-7.3f max|T - T_est| = %.1e\n', ...
             model.dt, columns(Tbe)-1, (columns(Tbe)-1)*model.dt, max(abs(Tbe(:,end)-Test)));
    end
    xl = xlim;
    plot(xl, Test(iC)*[1 1], 'k:');
    leyenda{end+1} = 'estacionario';
    hold off; grid on;
    xlabel('t'); ylabel(sprintf('T(%g, %g)', xc, yc));
    legend(leyenda, 'location', 'southeast');
    title('Transitorio en el centro');
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
