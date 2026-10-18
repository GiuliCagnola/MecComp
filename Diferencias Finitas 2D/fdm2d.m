function [PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model)
    % - xnode : matriz de pares (x,y). Cada fila contiene las coordenadas
    %   x e y de un nodo de la malla. El número de fila (índice) corresponde
    %   al número de nodo.
    % - icone : matriz de conectividad. Cada fila contiene los índices de
    %   los nodos que integran un elemento, típicamente los 4 nodos que forman
    %   un cuadrángulo.

    % Inicialización de variables principales del sistema
    [K,F] = fdm2d_initialize(model.nnodes);
    
    % Armado de la matriz de vecindad
    [neighb] = fdm2d_neighbors(icone);
    
    % Ensamble de coeficientes del sistema
    [K,F] = fdm2d_gen_system(K,F,xnode,neighb,model.k,model.c,model.G);

    if isfield(model, 'cb_orden') && model.cb_orden == 1
        % Neumann y Robin de 1er orden (sin nodo ficticio), sólo estacionario
        [K,F] = fdm2d_cb_orden1(K,F,xnode,neighb,model.k,NEU,ROB);
    else
        % Ensamble de nodos frontera Neumann
        [F] = fdm2d_neumann(F,xnode,neighb,NEU);

        % Ensamble de nodos frontera Robin
        [K,F] = fdm2d_robin(K,F,xnode,neighb,ROB);
    end
    
    % Ensamble de nodos frontera Dirichlet
    [K,F] = fdm2d_dirichlet(K,F,DIR);
    
    % Los esquemas temporales necesitan los nodos Dirichlet para fijar su
    % valor en cada paso de tiempo
    model.DIR = reshape(DIR, [], 2);    % 0x2 si no hay nodos Dirichlet

    % Resolución del sistema lineal de ecuaciones
    [PHI,Q] = fdm2d_solve(K,F,xnode,neighb,model);
end

