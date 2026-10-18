function [T, Thist, t] = elemFinitos(xnode, model, cb, et)
% elemFinitos  Transferencia de calor 1D por el Método de Elementos Finitos
%
%   rho*cp dT/dt = d/dx(k dT/dx) - c T + G      en [xnode(1), xnode(end)]
%
% ENTRADAS
%   xnode : vector de coordenadas nodales (ordenadas de izquierda a derecha)
%   model : struct con k, c, rho, cp, G
%             k, c y G pueden ser escalares o function_handle @(x) ...
%             model.orden (opcional): 1 = elementos lineales (por defecto)
%                                     2 = elementos cuadráticos: xnode debe
%                                     tener 2*Nelem+1 nodos; el elemento e usa
%                                     los nodos 2e-1, 2e, 2e+1
%   cb    : matriz 2x3, fila 1 = extremo izquierdo, fila 2 = extremo derecho
%             [1 T    -1]   Dirichlet
%             [2 q    -1]   Neumann, q = -k dT/dn (flujo SALIENTE)
%             [3 h  Tinf]   Robin,   q = h (T - Tinf)
%   et    : esquema temporal, escalar o vector [esquema dt tfinal]
%             esquema: 0 = estacionario, 1 = Forward Euler,
%                      2 = Backward Euler, 3 = Crank-Nicolson
%             dt     (opcional): por defecto 0.9*dt crítico de Forward Euler
%             tfinal (opcional): si no se da, se avanza hasta el estacionario
%
% SALIDAS
%   T     : temperaturas nodales (columna). En el transitorio, el último paso.
%   Thist : historia, una columna por paso de tiempo (transitorio)
%   t     : instantes de tiempo de cada columna de Thist
%
% Condición inicial del transitorio: T(x,0) = 0 (salvo nodos Dirichlet).

    if nargin < 4
        et = 0;
    end
    xnode = xnode(:);
    n = numel(xnode);
    orden = 1;
    if isfield(model, 'orden')
        orden = model.orden;
    end
    if orden == 2 && mod(n-1, 2) ~= 0
        error('Con elementos cuadráticos la cantidad de nodos debe ser impar.');
    end
    nen = orden + 1;                         % nodos por elemento
    ne  = (n-1)/orden;                       % cantidad de elementos
    icone = (1:orden:n-orden)' + (0:orden);  % conectividad ne x nen

    % coeficientes como funciones de x
    k     = comoFuncion(model.k);
    c     = comoFuncion(model.c);
    G     = comoFuncion(model.G);
    rhocp = model.rho*model.cp;

    % ── ensamble: K (difusión), M (reacción), C (capacidad), F (fuente) ──
    % integración de Gauss con 3 puntos (exacta hasta grado 5)
    xg = [-sqrt(3/5) 0 sqrt(3/5)];
    wg = [5 8 5]/9;
    K = zeros(n);  M = zeros(n);  C = zeros(n);  F = zeros(n,1);
    for e = 1 : ne
        nod = icone(e,:);
        xe = xnode(nod);
        Ke = zeros(nen);  Me = zeros(nen);  Ce = zeros(nen);  Fe = zeros(nen,1);
        for g = 1 : numel(xg)
            [N, dN] = funcionesForma(orden, xg(g));
            J   = dN*xe;                     % jacobiano dx/dxi
            x   = N*xe;                      % punto de Gauss en coordenadas reales
            dNx = dN/J;                      % derivadas respecto de x
            w   = wg(g)*J;
            Ke = Ke + w*k(x)*(dNx'*dNx);
            Me = Me + w*c(x)*(N'*N);
            Ce = Ce + w*rhocp*(N'*N);
            Fe = Fe + w*G(x)*N';
        end
        K(nod,nod) = K(nod,nod) + Ke;
        M(nod,nod) = M(nod,nod) + Me;
        C(nod,nod) = C(nod,nod) + Ce;
        F(nod)     = F(nod)     + Fe;
    end

    % ── condiciones de borde: A T = b, con A = K + M + H, b = F + Q ──
    A = K + M;
    b = F;
    extremos = [1 n];
    fijo = [];  valFijo = [];
    for s = 1 : 2
        P = extremos(s);
        switch cb(s,1)
            case 1      % Dirichlet: se impone después, sobre el sistema reducido
                fijo(end+1) = P;  valFijo(end+1) = cb(s,2);
            case 2      % Neumann: Q = -q
                b(P) = b(P) - cb(s,2);
            case 3      % Robin: H = h, Q = h*Tinf
                A(P,P) = A(P,P) + cb(s,2);
                b(P)   = b(P)   + cb(s,2)*cb(s,3);
            otherwise
                error('Tipo de condición de borde desconocido: %d', cb(s,1));
        end
    end
    libre = setdiff(1:n, fijo);
    valFijo = valFijo(:);

    % sistema reducido: la columna de los nodos fijos pasa al lado derecho
    Al = A(libre,libre);
    bl = b(libre) - A(libre,fijo)*valFijo;

    T = zeros(n,1);
    T(fijo) = valFijo;

    esquema = et(1);
    if esquema == 0
        T(libre) = Al \ bl;
        Thist = T;  t = 0;
        return
    end

    % ── transitorio: método theta ──
    if rhocp == 0
        error('rho*cp = 0: el problema no tiene término temporal. Usar et = 0.');
    end
    switch esquema
        case 1,  theta = 0;      % Forward Euler
        case 2,  theta = 1;      % Backward Euler
        case 3,  theta = 0.5;    % Crank-Nicolson
        otherwise, error('esquema debe ser 0, 1, 2 o 3');
    end
    Cl = C(libre,libre);

    % dt crítico de Forward Euler: 2 / lambda_max(C^-1 A)
    dtCrit = 2 / max(real(eig(Al, Cl)));
    if numel(et) >= 2 && et(2) > 0
        dt = et(2);
    else
        dt = 0.9*dtCrit;
    end
    if theta == 0 && dt > dtCrit
        warning('Forward Euler con dt = %g > dt crítico = %g: el esquema es inestable.', dt, dtCrit);
    end

    % (C + theta*dt*A) T^{n+1} = (C - (1-theta)*dt*A) T^n + dt*b
    Kth = Cl + theta*dt*Al;
    Rth = Cl - (1-theta)*dt*Al;
    [L, U, Pp] = lu(Kth);              % Kth no cambia: se factoriza una vez

    if numel(et) >= 3 && et(3) > 0
        nPasos = round(et(3)/dt);       % hasta tfinal
        hastaEstacionario = false;
    else
        nPasos = 1e6;                    % hasta el estacionario
        hastaEstacionario = true;
    end

    Tl = zeros(numel(libre),1);          % T(x,0) = 0
    Thist = zeros(n, min(nPasos,1000)+1);
    Thist(:,1) = T;
    for p = 1 : nPasos
        Tl1 = U \ (L \ (Pp*(Rth*Tl + dt*bl)));
        T(libre) = Tl1;
        if p+1 > columns(Thist)          % se agranda la historia de a bloques
            Thist(:, end+1:2*end) = 0;
        end
        Thist(:,p+1) = T;
        % estacionario: la velocidad de cambio es despreciable
        cambio = norm(Tl1 - Tl, inf)/dt;
        Tl = Tl1;
        if !all(isfinite(Tl))
            warning('La solución diverge (paso %d).', p);
            break;
        end
        if hastaEstacionario && cambio < 1e-6*max(1, norm(T, inf))
            break;
        end
    end
    Thist = Thist(:, 1:p+1);
    t = (0:p)*dt;
end

function f = comoFuncion(v)
    if isa(v, 'function_handle')
        f = v;
    else
        f = @(x) v;
    end
end

function [N, dN] = funcionesForma(orden, xi)
% Funciones de forma en el elemento de referencia xi en [-1, 1]
    if orden == 1
        N  = [(1-xi)/2, (1+xi)/2];
        dN = [-1/2, 1/2];
    else
        N  = [xi*(xi-1)/2, 1-xi^2, xi*(xi+1)/2];
        dN = [xi-1/2, -2*xi, xi+1/2];
    end
end
