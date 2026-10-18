function [PHI, Q] = fdm2d_implicit(K,F,xnode,neighb,model,dt)
% Descripción: módulo para resolver el sistema lineal de ecuaciones utilizando esquema
% temporal implícito (Backward Euler):
%   (rho*cp/dt*I + K) PHI^{n+1} = F + rho*cp/dt*PHI^n

% Entrada:
% * K: matriz del sistema (difusión + reacción)
% * F: vector de flujo térmico.
% * xnode: matriz de nodos con pares (x,y) representando las coordenadas de cada nodo
% de la malla.
% * neighb: matriz de vecindad.
% * model: struct con todos los datos del modelo (constantes, esquema numérico, etc.)
% * dt: paso temporal (arbitrario, el esquema es incondicionalmente estable).

% Salida:
% * PHI: matriz solución. Cada elemento del vector representa un valor escalar
%   asociado a cada nodo de la malla, y su posición dentro del vector depende de
%   cómo se especificó cada nodo en xnode. Cada columna representa una iteración
%   del esquema temporal (en total nit columnas).
% * Q: matriz de flujo de calor. Para cada nodo se halla un vector bidimensional
%   de flujo de calor, representado por un par (Qx,Qy). Cada par de columnas
%   representa una iteración del esquema temporal (en total 2×nit columnas).
% ----------------------------------------------------------------------
    printf('Paso de tiempo implícito: dt = %g\n', dt);
    nt = model.maxit;
    N = size(xnode,1);
    tol = model.tol;
    a = (model.rho*model.cp)/dt;

    % Nodos Dirichlet: fila identidad en Ks y valor impuesto en B
    nodosD = model.DIR(:,1);
    valD = model.DIR(:,2);

    Ks = a*speye(N) + K;
    Ks(nodosD,:) = 0;
    Ks(sub2ind([N N], nodosD, nodosD)) = 1;

    % Ks no cambia entre pasos: se factoriza una sola vez (LU dispersa)
    [L, U, Pr, Pc] = lu(Ks);

    PHIn = model.PHI_n;
    PHIn(nodosD) = valD;

    % Historia prealocada
    PHI = zeros(N, nt+1);
    Q = zeros(N, 2*(nt+1));
    PHI(:,1) = PHIn;
    Q(:,1:2) = fdm2d_flux(PHIn,neighb,xnode,model.k);

    convergio = false;
    for i=1:nt
        B = F + a*PHIn;
        B(nodosD) = valD;
        PHIa = Pc*(U\(L\(Pr*B)));

        err = norm(PHIa-PHIn,2)/norm(PHIa,2);
        PHIn = PHIa;
        PHI(:,i+1) = PHIa;
        Q(:,2*i+1:2*i+2) = fdm2d_flux(PHIa,neighb,xnode,model.k);
        if (err < tol)
            convergio = true;
            break;
        endif
    endfor

    % Descartamos las columnas no usadas
    PHI = full(PHI(:,1:i+1));
    Q = Q(:,1:2*(i+1));

    if convergio
        printf('Método implícito terminado por tolerancia de error (%d pasos).\n', i);
    else
        printf('Método implícito terminado por límite de iteraciones (%d pasos).\n', i);
    endif
endfunction
