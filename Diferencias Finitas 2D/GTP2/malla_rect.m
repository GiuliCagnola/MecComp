function [xnode, icone] = malla_rect(x1, y1, x2, y2, nx, ny, alfa)
% Malla estructurada de cuadrángulos sobre el rectángulo [x1,x2]x[y1,y2]
% con nx x ny intervalos.
% alfa = 0: malla uniforme.
% alfa > 0 (y < 1): malla no uniforme, más fina cerca de los bordes,
%   x = x1 + Lx*(s - alfa*sin(2*pi*s)/(2*pi)),  s uniforme en [0,1].
%   Al duplicar nx los nodos de la malla gruesa siguen estando en la fina.
% Numeración: primero en x (i), después en y (j). icone en sentido
% antihorario empezando por el vértice inferior izquierdo.
    if nargin < 7
        alfa = 0;
    end
    mapa = @(s) s - alfa*sin(2*pi*s)/(2*pi);
    xs = x1 + (x2-x1)*mapa(linspace(0,1,nx+1));
    ys = y1 + (y2-y1)*mapa(linspace(0,1,ny+1));

    [X, Y] = ndgrid(xs, ys);
    xnode = [X(:) Y(:)];

    id = @(i,j) (j-1)*(nx+1) + i;
    icone = zeros(nx*ny, 4);
    e = 0;
    for j = 1 : ny
        for i = 1 : nx
            e = e + 1;
            icone(e,:) = [id(i,j) id(i+1,j) id(i+1,j+1) id(i,j+1)];
        end
    end
end
