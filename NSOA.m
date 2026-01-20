%_________________________________________________________________________%
%                                                                         %
%  Neural Swarm Optimization Algorithm (NSOA) source codes version 1.0    %
%                                                                         %
%  Developed in:	MATLAB 9.5 (R2018b)                                   %
%                                                                         %
%  Programmer:		Jie Cai                                               %
%                                                                         %
%  Original paper:	Jie Cai, Jiayao Chen                                  %
%					A novel metaheuristic algorithm for solving           %
%                   global optimization and engineering problems.         %
%                                                                         %
%_________________________________________________________________________%
function [Fbest, Xbest, CNVG] = NSOA(N, T, lb, ub, dim, fobj)
%% ============ ONLY EXPLICIT PARAMETERS ============
alpha = 0.12;
rho   = 0.60;
rho   = max(0, min(1, rho));

%% ============ SCALE & BASIC PREP ============
if numel(lb) == 1, lb = repmat(lb, 1, dim); end
if numel(ub) == 1, ub = repmat(ub, 1, dim); end
lb = reshape(lb, 1, []);
ub = reshape(ub, 1, []);

L = mean(ub - lb);

%% ============ INITIALIZATION ============
X = initialization(N, dim, ub, lb);
V = zeros(N, dim);
Fitness = inf(N, 1);
for i = 1:N
    Fitness(i) = fobj(X(i,:));
end
[Fbest, bestIdx] = min(Fitness);
Xbest = X(bestIdx, :);
CNVG  = zeros(T, 1);

Heading = randn(N, dim);
Heading = Heading ./ (sqrt(sum(Heading.^2, 2)) + 1e-12);

%% ============ MAIN LOOP ============
for t = 1:T
    s = t / T;

    % step = alpha * L * (1-s)^kappa, with kappa=1.8 embedded
    step = alpha * L * (1 - s)^1.8;

    %% --- 1) Core pull ---
    core_k = max(3, round(0.15 * N));   % rule (embedded)
    [~, idxCore] = mink(Fitness, min(core_k, N));
    swarm_core = mean(X(idxCore, :), 1);

    toward_core = swarm_core - X;
    toward_core = toward_core ./ (sqrt(sum(toward_core.^2, 2)) + 1e-12);

    %% --- 2) Best-buddy beacon ---
    buddy_beacon = zeros(N, dim);
    for i = 1:N
        nbrs = setdiff(randi(N, 1, 5), i);
        if isempty(nbrs)
            others = setdiff(1:N, i);
            nbrs = others(randi(numel(others)));
        end
        nbrs = nbrs(1:min(3, numel(nbrs)));
        [~, bi] = min(Fitness(nbrs));
        best_nbr = nbrs(bi);
        vec = X(best_nbr, :) - X(i, :);
        buddy_beacon(i, :) = vec ./ (norm(vec) + 1e-12);
    end

    %% --- 3) No-bump repulsion (density-aware radius, embedded) ---
    no_bump_r = (L * N^(-1 / max(1, dim)));   % cr=1.0 embedded
    r2 = no_bump_r^2;

    no_bump = zeros(N, dim);
    for i = 1:N
        diffs = X - X(i, :);
        D2 = sum(diffs.^2, 2);
        D2(i) = inf;
        too_close = D2 < r2;
        if any(too_close)
            vecs  = diffs(too_close, :);
            dists = sqrt(D2(too_close));
            no_bump(i, :) = sum(vecs ./ (dists + 1e-12), 1);
        end
    end
    no_bump = no_bump ./ (sqrt(sum(no_bump.^2, 2)) + 1e-12);

    %% --- 4) Velocity update + consensus heading ---
    g  = 1 + s;
    wC = rho + (1 - rho) * s;
    wC = max(0, min(1, wC));
    wB = 1 - wC;

    V = 0.4 * V + g * (wB * buddy_beacon + wC * toward_core + 0.5 * no_bump);

    dir_unit = V ./ (sqrt(sum(V.^2, 2)) + 1e-12);
    tau = 0.3 + 0.6 * (rho * s + (1 - rho) * s^2);
    tau = max(0.3, min(0.9, tau));

    Hd = (1 - tau) * Heading + tau * dir_unit;
    Heading = Hd ./ (sqrt(sum(Hd.^2, 2)) + 1e-12);

    Vmags = sqrt(sum(V.^2, 2));
    V = Heading .* Vmags;

    X = X + step * V;

    %% --- 5) Rare Levy probing (all constants embedded; only rho drives intensity) ---
    march_order = norm(mean(Heading, 1)); 
    p_eff = 0.08 * (1 - rho) * (1 - march_order);
    p_eff = max(0, min(0.08, p_eff));
    num_feeler = round(N * p_eff);

    lambda = 0.02 * L * (1 - rho) * (1 - 0.5 * march_order);

    if num_feeler > 0 && lambda > 0
        feeler_idx = randperm(N, num_feeler);

        beta = 1.5;
        sigma_u = (gamma(1 + beta) / ...
                  (beta * gamma((1 + beta)/2) * 2^((beta - 1)/2)))^(1/beta);

        for ii = 1:num_feeler
            i = feeler_idx(ii);
            u = randn(1, dim) * sigma_u;
            v = randn(1, dim);
            levy_step = u ./ (abs(v).^(1/beta) + 1e-12);
            X(i, :) = X(i, :) + (lambda * levy_step) .* randn(1, dim);
        end
    end

    %% --- 6) Boundary reflection ---
    for d = 1:dim
        low = X(:, d) < lb(d);
        up  = X(:, d) > ub(d);
        X(low, d) = lb(d);  V(low, d) = -0.5 * V(low, d);
        X(up, d)  = ub(d);  V(up, d)  = -0.5 * V(up, d);
    end

    %% --- 7) Fitness evaluation ---
    for i = 1:N
        f = fobj(X(i, :));
        Fitness(i) = f;
        if f < Fbest
            Fbest = f;
            Xbest = X(i, :);
        end
    end
    CNVG(t) = Fbest;
end

CNVG(CNVG == 0) = Fbest;
end

%% ================== Initialization ==================
function X = initialization(N, dim, ub, lb)
    if numel(ub) == 1
        X = rand(N, dim) * (ub - lb) + lb;
    else
        X = zeros(N, dim);
        for i = 1:dim
            X(:, i) = rand(N, 1) * (ub(i) - lb(i)) + lb(i);
        end
    end
end
