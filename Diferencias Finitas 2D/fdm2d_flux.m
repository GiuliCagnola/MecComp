function [Q] = fdm2d_flux(PHI,neighb,xnode,k)
% Flujo de calor Q = -k*grad(PHI) en cada nodo.
% La derivada en cada dirección se calcula con 3 nodos alineados (derivada del
% polinomio interpolante de Lagrange), lo que da orden dx^2 también en mallas
% no uniformes y en los bordes:
%   - nodo interior: W, P, E (centrada)
%   - nodo de borde: P y los dos siguientes hacia adentro (descentrada)
% Si en esa dirección sólo hay 2 nodos se usa la diferencia de 2 puntos.
% El cálculo está vectorizado (se llama en cada paso de tiempo).
    Qx = -k(:) .* derivada(PHI, xnode(:,1), neighb(:,4), neighb(:,2));
    Qy = -k(:) .* derivada(PHI, xnode(:,2), neighb(:,1), neighb(:,3));
    Q = [Qx, Qy];
end

function d = derivada(PHI, x, ant, sig)
% Derivada de PHI respecto de la coordenada x en todos los nodos.
% ant/sig: vecino anterior/siguiente en esa dirección (-1 si no existe).
    n = numel(x);
    P = (1:n)';
    PHI = full(PHI(:));

    % ternas de nodos (i0, i1, i2) usadas en cada nodo
    i0 = ant;  i1 = P;  i2 = sig;              % nodo interior: centrada
    bIni = (ant == -1);                        % primer nodo de la línea
    bFin = (sig == -1);                        % último nodo de la línea
    i0(bIni) = P(bIni);    i1(bIni) = sig(bIni);
    i2(bIni) = sig(max(sig(bIni),1));          % vecino del vecino
    i2(bIni & sig == -1) = -1;
    i0(bFin) = ant(max(ant(bFin),1));          % vecino del vecino
    i0(bFin & ant == -1) = -1;
    i1(bFin) = ant(bFin);  i2(bFin) = P(bFin);
    % (i0 de bFin queda en -1 si la línea tiene sólo 2 nodos)

    d = zeros(n,1);
    xp = x;

    % --- 3 puntos: derivada de las bases de Lagrange en xp ---
    t = (i0 ~= -1) & (i1 ~= -1) & (i2 ~= -1);
    x0 = x(i0(t)); x1 = x(i1(t)); x2 = x(i2(t)); xq = xp(t);
    w0 = ((xq-x1) + (xq-x2)) ./ ((x0-x1).*(x0-x2));
    w1 = ((xq-x0) + (xq-x2)) ./ ((x1-x0).*(x1-x2));
    w2 = ((xq-x0) + (xq-x1)) ./ ((x2-x0).*(x2-x1));
    d(t) = w0.*PHI(i0(t)) + w1.*PHI(i1(t)) + w2.*PHI(i2(t));

    % --- 2 puntos (línea con sólo 2 nodos): pendiente ---
    % falta i2 -> se usan (i0,i1); falta i0 -> se usan (i1,i2)
    a = i0;  b = i1;
    a(i0 == -1) = i1(i0 == -1);
    b(i0 == -1) = i2(i0 == -1);
    t2 = !t & (a ~= -1) & (b ~= -1);
    d(t2) = (PHI(b(t2)) - PHI(a(t2))) ./ (x(b(t2)) - x(a(t2)));
end
