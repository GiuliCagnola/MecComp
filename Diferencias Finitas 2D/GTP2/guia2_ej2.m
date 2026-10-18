close all; clear all; more off;
% =========================================================================
% GUIA DE TRABAJOS PRACTICOS N 2 - Diferencias Finitas 2D - Ejercicio 2
%
% Dominio en L: (0,0) (3,0) (3,2) (1,2) (1,6) (0,6).  k = 1, c = 0.
%   x = 0        : T = 20
%   y = 6        : T = 0
%   x = 1 (y>2)  : q = 0          y = 2 (x>1) : q = 0
%   x = 3        : q = 0          y = 0       : Robin h = 2, Tinf = 50
%
% a) Estacionario con G = 100.
% b) Transitorio con rho*cp = 1, dt = 0.05 (Backward Euler), tol = 1e-5.
% c) Estacionario con G = 100 solo en 0<=x<=1, 0<=y<=3 (G = 0 en el resto).
% d) Item b, pero la fuente se anula a los 2 segundos.
% Cada item se resuelve con distintos refinamientos (h = 0.5 ... 0.0625).
%
% Esquinas: en (0,6) se juntan T=20 y T=0 (se toma el promedio; ese nodo no
% interviene en el stencil de ningun nodo interior). En (0,0) y (1,6) manda
% el Dirichlet.
% =========================================================================
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

hRef  = [0.5 0.25 0.125 0.0625];   % refinamientos
hGraf = 0.125;                     % malla para graficar
dt    = 0.05;  tol = 1e-5;
puntos = [0.5 4; 2 1; 0.5 1];      % puntos de control
nomPts = {'(0.5, 4)', '(2, 1)', '(0.5, 1)'};

function [DIR, NEU, ROB] = cbL(xnode, icone)
  s(1) = struct('dir', 4, 'pos', 0, 'rango', [0 6], 'cb', [1 20 0]);    % x = 0
  s(2) = struct('dir', 3, 'pos', 6, 'rango', [0 1], 'cb', [1 0 0]);     % y = 6
  s(3) = struct('dir', 2, 'pos', 1, 'rango', [2 6], 'cb', [2 0 0]);     % x = 1
  s(4) = struct('dir', 3, 'pos', 2, 'rango', [1 3], 'cb', [2 0 0]);     % y = 2
  s(5) = struct('dir', 2, 'pos', 3, 'rango', [0 2], 'cb', [2 0 0]);     % x = 3
  s(6) = struct('dir', 1, 'pos', 0, 'rango', [0 3], 'cb', [3 2 50]);    % y = 0
  [DIR, NEU, ROB] = cb_segmentos(xnode, icone, s);
endfunction

function [xnode, icone, DIR, NEU, ROB, model] = armarL(h, fuente)
  [xnode, icone] = malla_L(h);
  [DIR, NEU, ROB] = cbL(xnode, icone);
  nn = rows(xnode);
  model.nnodes = nn;
  model.k = ones(nn,1);  model.c = zeros(nn,1);
  model.G = fuente(xnode(:,1), xnode(:,2));
  model.rho = 1;  model.cp = 1;
  model.ts = 2;  model.PHI_n = zeros(nn,1);
endfunction

function idx = buscar(xnode, pts)
  idx = zeros(rows(pts),1);
  for i = 1 : rows(pts)
    [~, idx(i)] = min((xnode(:,1)-pts(i,1)).^2 + (xnode(:,2)-pts(i,2)).^2);
  end
endfunction

% Estacionario en todas las mallas + orden observado (mallas sucesivas
% comparadas en los nodos de la malla mas gruesa)
function [Tgraf, xg, ig, Tpts] = refinamiento(hRef, hGraf, fuente, puntos, nomPts, titulo)
  printf('\n%s\n', titulo);
  printf('%8s %7s', 'h', 'nodos');  printf('   T%-10s', nomPts{:});  printf('\n');
  Tpts = zeros(numel(hRef), rows(puntos));
  for j = 1 : numel(hRef)
    [xnode, icone, DIR, NEU, ROB, model] = armarL(hRef(j), fuente);
    T = resolver(xnode, icone, DIR, NEU, ROB, model);
    Tpts(j,:) = T(buscar(xnode, puntos));
    if j == 1
      xg0 = xnode;                    % nodos de la malla mas gruesa
    end
    sub{j} = T(buscar(xnode, xg0));
    if hRef(j) == hGraf
      Tgraf = T; xg = xnode; ig = icone;
    end
    printf('%8.4f %7d', hRef(j), rows(xnode));  printf('   %-11.4f', Tpts(j,:));  printf('\n');
  end
  for j = 1 : numel(hRef)-1
    dif(j) = sqrt(mean((sub{j} - sub{j+1}).^2));
  end
  p = log2(dif(1:end-1)./dif(2:end));
  printf('RMS(T_h - T_h/2): %s\n', num2str(dif, '%10.3e'));
  printf('orden observado p: %s\n', num2str(p, '%8.2f'));

  loglog(hRef(1:end-1), dif, 'bo-', 'linewidth', 1.5); hold on;
  hh = hRef(1:end-1);
  loglog(hh, dif(end)*(hh/hh(end)).^2, 'k--');
  hold off; grid on;
  xlabel('h'); ylabel('RMS(T_h - T_{h/2})');
  legend(sprintf('p = %.2f', p(end)), 'pendiente 2', 'location', 'northwest');
  title('Convergencia');
endfunction

G100   = @(x,y) 100*ones(size(x));
Gzona  = @(x,y) 100*((x <= 1+1e-10) & (y <= 3+1e-10));
G0     = @(x,y) zeros(size(x));

% ======================= a) Estacionario, G = 100 =======================
figure('Name', 'Ejercicio 2.a', 'Position', [50 50 1000 500]);
subplot(1,2,2);
[Ta, xg, ig] = refinamiento(hRef, hGraf, G100, puntos, nomPts, '=== a) Estacionario, G = 100 ===');
subplot(1,2,1);
graficar_campo(xg, ig, Ta, sprintf('2.a: estacionario, G = 100 (h = %g)', hGraf));

% ======================= b) Transitorio, Backward Euler =======================
printf('\n=== b) Transitorio, rho*cp = 1, dt = %g, tol = %g (Backward Euler) ===\n', dt, tol);
printf('%8s %7s %7s %9s %16s\n', 'h', 'nodos', 'pasos', 't final', 'max|T - T_est|');
figure('Name', 'Ejercicio 2.b', 'Position', [50 50 1000 500]);
subplot(1,2,2); hold on;
estilos = {'-', '--', ':'};  colores = {'b', 'r'};  leyenda = {};   % estilo: malla, color: punto
for j = 2 : numel(hRef)
  h = hRef(j);
  [xnode, icone, DIR, NEU, ROB, model] = armarL(h, G100);
  Test = resolver(xnode, icone, DIR, NEU, ROB, model);
  model.ts = 1;  model.dt = dt;  model.maxit = 10000;  model.tol = tol;
  Tb = resolver(xnode, icone, DIR, NEU, ROB, model);
  np = columns(Tb) - 1;
  printf('%8.4f %7d %7d %9.2f %16.2e\n', h, rows(xnode), np, np*dt, max(abs(Tb(:,end)-Test)));
  ip = buscar(xnode, puntos);
  t = (0:np)*dt;
  for m = 1 : 2
    plot(t, Tb(ip(m),:), [colores{m} estilos{j-1}], 'linewidth', 1.5);
    leyenda{end+1} = sprintf('T%s, h = %g', nomPts{m}, h);
  end
  if h == hGraf
    TbGraf = Tb; xb = xnode; ib = icone; tb = t; ipb = ip;
  end
end
hold off; grid on; xlabel('t'); ylabel('T');
legend(leyenda, 'location', 'east'); title('2.b: evolución en los puntos de control');
subplot(1,2,1);
graficar_campo(xb, ib, TbGraf(:,end), sprintf('2.b: estado final, t = %.2f (h = %g)', tb(end), hGraf));

% ======================= c) Estacionario, fuente localizada =======================
figure('Name', 'Ejercicio 2.c', 'Position', [50 50 1000 500]);
subplot(1,2,2);
[Tc, xg, ig] = refinamiento(hRef, hGraf, Gzona, puntos, nomPts, ...
                            '=== c) Estacionario, G = 100 en 0<=x<=1, 0<=y<=3 ===');
subplot(1,2,1);
graficar_campo(xg, ig, Tc, sprintf('2.c: G = 100 en x<=1, y<=3 (h = %g)', hGraf));
hold on; plot([0 1 1 0 0], [0 0 3 3 0], 'w--', 'linewidth', 1.5); hold off;

% ======================= d) Fuente apagada a los 2 s =======================
printf('\n=== d) Item b con la fuente anulada en t = 2 ===\n');
printf('%8s %7s %7s %9s %12s %16s\n', 'h', 'nodos', 'pasos', 't final', 'T max(t=2)', 'max|T - T_est|');
nOn = round(2/dt);                 % pasos con fuente encendida
figure('Name', 'Ejercicio 2.d', 'Position', [50 50 1000 500]);
subplot(1,2,2); hold on; leyenda = {};
for j = 2 : numel(hRef)
  h = hRef(j);
  % etapa 1: 0 < t <= 2 con G = 100 (tol = 0 -> hace exactamente nOn pasos)
  [xnode, icone, DIR, NEU, ROB, model] = armarL(h, G100);
  model.ts = 1;  model.dt = dt;  model.maxit = nOn;  model.tol = 0;
  T1 = resolver(xnode, icone, DIR, NEU, ROB, model);
  % etapa 2: t > 2 con G = 0, partiendo del estado en t = 2
  model.G = G0(xnode(:,1), xnode(:,2));
  model.PHI_n = T1(:,end);  model.maxit = 10000;  model.tol = tol;
  T2 = resolver(xnode, icone, DIR, NEU, ROB, model);
  Td = [T1 T2(:,2:end)];
  % estacionario con G = 0 (al que deberia llegar)
  model.ts = 2;
  Test0 = resolver(xnode, icone, DIR, NEU, ROB, model);
  np = columns(Td) - 1;
  printf('%8.4f %7d %7d %9.2f %12.4f %16.2e\n', h, rows(xnode), np, np*dt, max(T1(:,end)), max(abs(Td(:,end)-Test0)));
  ip = buscar(xnode, puntos);
  t = (0:np)*dt;
  for m = 1 : 2
    plot(t, Td(ip(m),:), [colores{m} estilos{j-1}], 'linewidth', 1.5);
    leyenda{end+1} = sprintf('T%s, h = %g', nomPts{m}, h);
  end
  if h == hGraf
    TdGraf = Td; xd = xnode; id = icone; td = t; Test0Graf = Test0;
  end
end
yl = ylim; plot([2 2], yl, 'k:'); leyenda{end+1} = 'se apaga G';
hold off; grid on; xlabel('t'); ylabel('T');
legend(leyenda, 'location', 'northeast'); title('2.d: evolución en los puntos de control');
subplot(1,2,1);
[~, i2] = min(abs(td - 2));
graficar_campo(xd, id, TdGraf(:,i2), sprintf('2.d: campo en t = 2, al apagar G (h = %g)', hGraf));

% Si se corre fuera de la GUI (p.ej. "octave archivo.m" desde VS Code),
% Octave cierra las figuras al terminar el script: esperamos a que el
% usuario las cierre.
if !isguirunning()
  disp('Cerrar las figuras para terminar...');
  for h = findobj('type', 'figure')'
    waitfor(h);
  endfor
endif
