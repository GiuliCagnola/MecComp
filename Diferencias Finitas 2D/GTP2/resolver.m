function [PHI, Q, dt] = resolver(xnode, icone, DIR, NEU, ROB, model)
% Llama a fdm2d sin mostrar sus mensajes. Devuelve también el paso de tiempo
% usado (NaN en estado estacionario).
    texto = evalc('[PHI, Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);');
    tok = regexp(texto, 'dt = ([-+0-9.eE]+)', 'tokens', 'once');
    if isempty(tok)
        dt = NaN;
    else
        dt = str2double(tok{1});
    end
end
