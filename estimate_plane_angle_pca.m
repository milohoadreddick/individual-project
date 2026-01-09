function [tc, theta] = estimate_plane_angle_pca(t, x, y, winN, stepN)

N = numel(t);
starts = 1:stepN:(N - winN + 1);
M = numel(starts);

tc = zeros(M,1);
theta = zeros(M,1);

v_prev = []; 

for i = 1:M
    s = starts(i);
    e = s + winN - 1;

    Xw = x(s:e) - mean(x(s:e));
    Yw = y(s:e) - mean(y(s:e));

    C = [mean(Xw.^2), mean(Xw.*Yw);
         mean(Xw.*Yw), mean(Yw.^2)];

    [V,D] = eig(C);
    [~,imax] = max(diag(D));
    v = V(:,imax);

    if ~isempty(v_prev)
        if dot(v, v_prev) < 0
            v = -v;
        end
    end
    v_prev = v;

    theta(i) = atan2(v(2), v(1));
    tc(i) = t(round((s+e)/2));
end
end