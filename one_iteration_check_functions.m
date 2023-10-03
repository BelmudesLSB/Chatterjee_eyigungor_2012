function c = consumption(q, b_grid, current_y, current_b, lambda, z)

    nb = length(b_grid);
    aux = zeros(nb);

    for x=1:nb
        aux(x) = -q(current_y * nb + x) * (b_grid(x) - (1-lambda) * b_grid(current_b)) + (lambda + z * (1-lambda)) * b_grid(current_b);
    end

    c = aux;

end