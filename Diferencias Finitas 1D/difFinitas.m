function [T] = difFinitas(xnode, model, cb, et, graficar)
  % graficar (opcional, por defecto true): con false no dibuja ni imprime
  if nargin < 5
    graficar = true;
  endif

  % inicializacion de variables %
  xnode = xnode(:);
  N = length(xnode);
  K = zeros(N,N);
  F = zeros(N,1);

  % parametros del modelo %
  rho_cp = model.rho*model.cp;
  k = model.k;
  cr = model.c;

  % ecuaciones para nodos interiores %
  for i=2:N-1
    dx_izq = xnode(i)-xnode(i-1);
    dx_der = xnode(i+1)-xnode(i);

    a = 2/(dx_izq*(dx_der+dx_izq));
    b = 2/(dx_der*dx_izq);
    c = 2/(dx_der*(dx_der+dx_izq));

    A = (k*a);
    B = -(k*b+cr);
    C = (k*c);

    K(i,i-1:i+1) = [A B C];
    F(i) = model.G(i);
  endfor

  % Cond. de borde en x=L1 %
  dx = xnode(2)-xnode(1);
  d2x = dx*dx;
  if cb(1,1) == 2
    # Cond neumann
    q = cb(1,2);
    B = -((2*k)/d2x+cr);
    C = (2*k)/d2x;
    K(1,1:2) = [B C];
    F(1) = (model.G(1)-2*q/dx);
  elseif cb(1,1) == 3
    # Cond robin
    h = cb(1,2);
    Tinf = cb(1,3);
    B = -((2*h)/dx+(2*k)/d2x+cr);
    C = (2*k)/d2x;
    K(1,1:2) = [B C];
    F(1) = (model.G(1)+(2*h*Tinf)/dx);
  endif

  % Cond. de borde en x=L2 %
  dx = xnode(N)-xnode(N-1);
  d2x = dx*dx;
  if cb(2,1) == 2
    # Cond neumann
    q = cb(2,2);
    A = (2*k)/d2x;
    B = -((2*k)/d2x+cr);
    K(N,N-1:N) = [A B];
    F(N) = (model.G(N)-2*q/dx);
  elseif cb(2,1) == 3
    # Cond robin
    h = cb(2,2);
    Tinf = cb(2,3);
    A = (2*k)/d2x;
    B = -((2*h)/dx+(2*k)/d2x+cr);
    K(N,N-1:N) = [A B];
    F(N) = (model.G(N)+(2*h*Tinf)/dx);
  endif

  % nodos con cond. dirichlet y sus valores %
  nodosD = [];
  valD = [];
  if cb(1,1) == 1
    nodosD(end+1) = 1;  valD(end+1) = cb(1,2);
  endif
  if cb(2,1) == 1
    nodosD(end+1) = N;  valD(end+1) = cb(2,2);
  endif

  if et(1) == 0 # estado estacionario
    F *= -1;
    [K, F] = dirchlet(K, F, cb);
    T = K\F;
    if graficar
      plot(xnode, T, 'b-x')
      disp('Estado estacionario');
    endif
    return
  endif

  % ----- esquemas transitorios ----- %
  maxIt = et(2); tol = et(3);

  if et(1) == 1 # esquema explicito
    # Estabilidad de Euler explicito: dt*|lambda_max|/rho_cp <= 2.
    # Por Gershgorin |lambda_max| <= 2*max|K_ii|, entonces
    # dt <= rho_cp/max|K_ii|. Esto considera la malla no uniforme
    # (el dx mas chico), el termino reactivo y los bordes Robin/Neumann.
    # En malla uniforme y c=0 se reduce a dx^2/(2*alpha).
    dt = 0.9*rho_cp/max(abs(diag(K))); # 0.9*dt_critico
  else # esquema implicito
    dt = et(4);
    Kimp = (rho_cp/dt)*eye(N,N) - K;
    # fila identidad en los nodos dirichlet
    Kimp(nodosD,:) = 0;
    Kimp(sub2ind([N N], nodosD, nodosD)) = 1;
    # Kimp no cambia entre pasos: se factoriza una sola vez
    [Lf, Uf, Pf] = lu(Kimp);
  endif
  if graficar
    printf('dt = %g\n', dt);
  endif

  # condicion inicial nula, salvo en los nodos dirichlet
  Tant = zeros(N,1);
  Tant(nodosD) = valD;

  T = zeros(N, maxIt+1); # historia de la solucion (prealocada)
  T(:,1) = Tant;

  # se actualiza la figura ~50 veces en total (no en cada iteracion)
  nPlot = max(1, floor(maxIt/50));
  if graficar
    hLinea = plot(xnode, Tant, '*-b');
  endif

  convergio = false;
  for i=1:maxIt
    if et(1) == 1
      Tact = (dt/rho_cp)*(K*Tant+F)+Tant;
    else
      Fimp = F + (rho_cp/dt)*Tant;
      Fimp(nodosD) = valD;
      Tact = Uf\(Lf\(Pf*Fimp));
    endif

    # fijamos la/las condiciones dirichlet
    Tact(nodosD) = valD;

    T(:,i+1) = Tact;
    err = norm(Tact-Tant,2)/norm(Tact,2);
    convergio = (err < tol);

    if graficar && (convergio || mod(i, nPlot) == 0 || i == maxIt)
      set(hLinea, 'YData', Tact);
      title(sprintf('nit: %d - t: %g - error: %e', i, i*dt, err));
      drawnow;
    endif

    Tant = Tact;
    if convergio
      break
    endif
  endfor

  T = T(:,1:i+1); # descartamos las columnas no usadas
  if !graficar
    return
  endif
  if convergio
    disp('Se alcanza el estado estacionario segun la tolerancia definida');
  else
    disp('Se alcanza el maximo de iteraciones sin llegar a un estado estacionario');
  endif

endfunction
