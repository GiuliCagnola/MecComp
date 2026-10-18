function [PHI_vec, Q_vec] = fdm2d_explicit(K,F,xnode,neighb,model,dt)
% Descripción: módulo para resolver el problema transitorio con esquema
% explícito (Forward Euler):
%   rho*cp*(PHI^{n+1} - PHI^n)/dt = F - K*PHI^n
%
% Salida:
% * PHI_vec: una columna por paso de tiempo (la primera es la condición inicial).
% * Q_vec: un par de columnas (Qx,Qy) por paso de tiempo.
% ----------------------------------------------------------------------
    nnodes = model.nnodes;
    maxit = model.maxit;
    A = dt/(model.rho*model.cp);
    printf('Paso de tiempo explícito: dt = %g\n', dt);

    % Nodos Dirichlet: se fijan en la condición inicial y en cada paso
    nodosD = model.DIR(:,1);
    valD = model.DIR(:,2);

    % Matriz de iteración (no cambia entre pasos)
    M = speye(nnodes) - A*K;
    AF = A*F;

    PHI_n = model.PHI_n;
    PHI_n(nodosD) = valD;

    % Historia prealocada
    PHI_vec = zeros(nnodes, maxit+1);
    Q_vec = zeros(nnodes, 2*(maxit+1));
    PHI_vec(:,1) = PHI_n;
    Q_vec(:,1:2) = fdm2d_flux(PHI_n,neighb,xnode,model.k);

    convergio = false;
    for n = 1 : maxit
        PHI = AF + M*PHI_n;
        PHI(nodosD) = valD;

        err = norm(PHI-PHI_n,2)/norm(PHI,2);

        PHI_vec(:,n+1) = PHI;
        Q_vec(:,2*n+1:2*n+2) = fdm2d_flux(PHI,neighb,xnode,model.k);
        PHI_n = PHI;

        if err < model.tol
            convergio = true;
            break;
        end
    end

    % Descartamos las columnas no usadas
    PHI_vec = full(PHI_vec(:,1:n+1));
    Q_vec = Q_vec(:,1:2*(n+1));

    if convergio
        printf('Método explícito terminado por tolerancia de error (%d pasos).\n', n);
    else
        printf('Método explícito terminado por límite de iteraciones (%d pasos).\n', n);
    end
end
