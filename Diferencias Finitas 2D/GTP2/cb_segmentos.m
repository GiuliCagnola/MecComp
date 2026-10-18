function [DIR, NEU, ROB] = cb_segmentos(xnode, icone, segs)
% Arma DIR, NEU y ROB a partir de una lista de segmentos de borde.
% Cada segmento es un struct con:
%   dir   : normal exterior del borde: 1 (S), 2 (E), 3 (N), 4 (W)
%   pos   : coordenada fija del borde (y para S/N, x para E/W)
%   rango : [a b] intervalo de la otra coordenada
%   cb    : [tipo v1 v2]  tipo 1 Dirichlet (v1 = T)
%                         tipo 2 Neumann   (v1 = q saliente)
%                         tipo 3 Robin     (v1 = h, v2 = phi_inf)
% Un nodo recibe la condición sólo si no tiene vecino en la dirección de la
% normal (es decir, si es un nodo de borde para el stencil).
% Esquinas: si alguno de los bordes es Dirichlet manda el Dirichlet (si hay
% dos valores distintos se promedian); si no, se aplican las dos condiciones,
% cada una en su dirección.
    tol = 1e-10;
    neighb = fdm2d_neighbors(icone);
    nn = rows(xnode);
    valD = cell(nn,1);          % valores Dirichlet de cada nodo
    otras = cell(nn,1);         % [dir tipo v1 v2] Neumann/Robin de cada nodo

    for s = 1 : numel(segs)
        sg = segs(s);
        if any(sg.dir == [1 3])         % borde horizontal
            fija = xnode(:,2);  libre = xnode(:,1);
        else                            % borde vertical
            fija = xnode(:,1);  libre = xnode(:,2);
        end
        nodos = find(abs(fija - sg.pos) < tol & ...
                     libre >= sg.rango(1) - tol & libre <= sg.rango(2) + tol);
        for P = nodos'
            if neighb(P, sg.dir) ~= -1
                continue;               % tiene vecino hacia afuera: no es borde
            end
            if sg.cb(1) == 1
                valD{P}(end+1) = sg.cb(2);
            else
                otras{P}(end+1,:) = [sg.dir sg.cb(1:3)];
            end
        end
    end

    DIR = zeros(0,2);  NEU = zeros(0,3);  ROB = zeros(0,4);
    for P = 1 : nn
        if !isempty(valD{P})
            DIR(end+1,:) = [P mean(valD{P})];
        else
            for r = 1 : rows(otras{P})
                o = otras{P}(r,:);
                if o(2) == 2
                    NEU(end+1,:) = [P o(3) o(1)];
                else
                    ROB(end+1,:) = [P o(3) o(4) o(1)];
                end
            end
        end
    end
end
