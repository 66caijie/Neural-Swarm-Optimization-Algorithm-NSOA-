function X = Initialization(N, dim, ub, lb)
    if numel(ub) == 1
        X = rand(N, dim) * (ub - lb) + lb;
    else
        X = zeros(N, dim);
        for i = 1:dim
            X(:, i) = rand(N, 1) * (ub(i) - lb(i)) + lb(i);
        end
    end
end

