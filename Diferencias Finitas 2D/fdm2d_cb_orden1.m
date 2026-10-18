function [K,F] = fdm2d_cb_orden1(K,F,xnode,neighb,k,NEU,ROB)
% Descripción: condiciones de borde Neumann y Robin con aproximación de
% 1er orden (diferencia de 2 puntos, sin nodo ficticio). La fila del nodo de
% borde se reemplaza por la ecuación de la condición de borde:
%   Neumann: -k dphi/dn = q             ->  k/d*(phiP - phiI) = -q
%   Robin:   -k dphi/dn = h(phi - phinf) ->  (k/d + h)*phiP - k/d*phiI = h*phinf
% donde I es el vecino interior en la dirección normal y d la distancia P-I.
% En una esquina con dos condiciones se suman ambas ecuaciones.
%
% NEU y ROB tienen el mismo formato que en fdm2d_neumann y fdm2d_robin.
% ----------------------------------------------------------------------
    NEU = reshape(NEU, [], 3);
    ROB = reshape(ROB, [], 4);

    % [nodo direccion tipo valor1 valor2]
    lista = [NEU(:,1) NEU(:,3) 2*ones(rows(NEU),1) NEU(:,2) zeros(rows(NEU),1);
             ROB(:,1) ROB(:,4) 3*ones(rows(ROB),1) ROB(:,2) ROB(:,3)];

    % vecino interior: borde S(1) -> N(3), E(2) -> W(4), N(3) -> S(1), W(4) -> E(2)
    opuesto = [3 4 1 2];
    usado = false(size(K,1),1);

    for i = 1 : rows(lista)
        P = lista(i,1);
        I = neighb(P, opuesto(lista(i,2)));
        d = norm(xnode(I,:) - xnode(P,:));

        if !usado(P)
            % se descarta la ecuación de la EDO de este nodo
            K(P,:) = 0;
            F(P) = 0;
            usado(P) = true;
        end

        K(P,P) = K(P,P) + k(P)/d;
        K(P,I) = K(P,I) - k(P)/d;
        if lista(i,3) == 2          % Neumann
            F(P) = F(P) - lista(i,4);
        else                        % Robin
            K(P,P) = K(P,P) + lista(i,4);
            F(P) = F(P) + lista(i,4)*lista(i,5);
        end
    end
end
