function graficar_campo(xnode, icone, T, titulo)
% Dibuja el campo escalar T sobre la malla (en el subplot/figura actual).
    M = rows(icone);
    tri = [icone(:,[1 2 3]); icone(:,[1 3 4])];
    patch('Vertices', xnode, 'Faces', tri, 'FaceVertexCData', T(:), ...
          'FaceColor', 'interp', 'EdgeColor', 'none');
    if M <= 400
        hold on;
        patch('Vertices', xnode, 'Faces', icone, 'FaceColor', 'none', ...
              'EdgeColor', [0.3 0.3 0.3]);
        hold off;
    end
    axis equal tight; colorbar; box on;
    xlabel('x'); ylabel('y');
    title(titulo);
end
