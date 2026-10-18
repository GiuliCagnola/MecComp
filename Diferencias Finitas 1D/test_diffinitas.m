clear; clc; close all;

# ---------- Malla ----------
L = 1;           # longitud de la barra [m]
N = 21;          # numero de nodos
xnode = linspace(0, L, N)';

# ---------- Modelo ----------
model.k   = 50;          # difusividad / conductividad [W/m K]
model.c   = 0;           # termino reactivo (sin reaccion)
model.rho = 7800;        # densidad [kg/m3] (acero, referencia)
model.cp  = 460;         # calor especifico [J/kg K]
model.G   = zeros(N,1);  # sin fuente distribuida

# ---------- Condiciones de borde ----------
# fila 1 = izquierdo, fila 2 = derecho
# columna 1: 1=Dirichlet, 2=Neumann, 3=Robin
cb = [1 100 -1;   # izquierdo: Dirichlet, T = 100 C
      1  20 -1];  # derecho:  Dirichlet, T = 20 C

# ========================================================
# CASO 1: Estado estacionario
# ========================================================
et = [0];
figure;
T_steady = difFinitas(xnode, model, cb, et);
disp("Solucion estacionaria (deberia ser una linea recta 100 -> 20):");
disp(T_steady');

# ========================================================
# CASO 2: Transitorio - Backward Euler
# ========================================================
# et(2)=max_it, et(3)=tol, et(4)=dt
# alpha = k/(rho*cp) ~ 1.4e-5 m2/s -> tiempo caracteristico L^2/(pi^2*alpha) ~ 7300 s,
# por eso dt tiene que ser del orden de decenas/cientos de segundos
et = [2, 2000, 1e-8, 100];
T_trans = difFinitas(xnode, model, cb, et);

figure;
hold on;
step = max(1, floor(size(T_trans,2)/8));  # graficar ~8 curvas intermedias
for j = 1:step:size(T_trans,2)
  plot(xnode, T_trans(:,j));
endfor
plot(xnode, T_steady, "k--x", "linewidth", 2); # referencia estacionaria
xlabel("x [m]"); ylabel("T [C]");
title("Evolucion temporal (Backward Euler) hacia el estacionario");
legend("evolucion...", "", "", "", "", "", "", "", "estacionario", "location", "eastoutside");
hold off;

printf("Ultimo estado transitorio vs estacionario (deberian casi coincidir):\n");
disp([T_trans(:,end), T_steady]);


# Si se corre fuera de la GUI (p.ej. "octave archivo.m" desde VS Code),
# Octave cierra las figuras al terminar el script: esperamos a que el
# usuario las cierre.
if !isguirunning()
  disp('Cerrar las figuras para terminar...');
  for h = findobj('type', 'figure')'
    waitfor(h);
  endfor
endif
