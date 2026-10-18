function [dt] = fdm2d_explicit_delta_t(K,model)
% Paso de tiempo para Forward Euler: PHI^{n+1} = PHI^n + dt/(rho*cp)*(F - K*PHI^n).
% Es estable si dt/(rho*cp)*lambda_max(K) <= 2. Por el teorema de Gershgorin,
% lambda_max(K) <= max_i sum_j |K(i,j)|, entonces
%   dt <= 2*rho*cp / max_i sum_j |K(i,j)|
% Este criterio considera automaticamente el dx/dy mas chico, el k mas grande
% (k variable), el termino reactivo y los bordes Robin.
% En malla uniforme (dx=dy), con k constante y sin reaccion se reduce a
% dt <= dx^2/(4*alpha), el criterio clasico en 2D.
% Las filas Dirichlet no cuentan: esos nodos se fijan en cada paso.

    filas = true(size(K,1),1);
    filas(model.DIR(:,1)) = false;

    radio = max(full(sum(abs(K(filas,:)),2)));
    dt = 0.9 * 2*model.rho*model.cp / radio;    % 0.9*dt_critico
end
