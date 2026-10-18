function [xnode, icone] = malla_L(h)
% Malla uniforme (paso h) del dominio en L del Ejercicio 2:
% vértices (0,0) (3,0) (3,2) (1,2) (1,6) (0,6).  1/h debe ser entero.
    tol = 1e-10;
    xs = 0:h:3;  ys = 0:h:6;
    nx = numel(xs);  ny = numel(ys);
    dentro = @(x,y) (x <= 1+tol) | (y <= 2+tol);

    % numeración de los nodos que quedan dentro del dominio
    num = zeros(nx, ny);
    xnode = [];
    for j = 1 : ny
        for i = 1 : nx
            if dentro(xs(i), ys(j))
                xnode(end+1,:) = [xs(i) ys(j)];
                num(i,j) = rows(xnode);
            end
        end
    end

    % celdas cuyo centro está dentro de la L
    icone = [];
    for j = 1 : ny-1
        for i = 1 : nx-1
            xc = (xs(i)+xs(i+1))/2;  yc = (ys(j)+ys(j+1))/2;
            if (xc < 1) || (yc < 2)
                icone(end+1,:) = [num(i,j) num(i+1,j) num(i+1,j+1) num(i,j+1)];
            end
        end
    end
end
