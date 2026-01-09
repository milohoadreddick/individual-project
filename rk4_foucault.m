function Z = rk4_foucault(t, z0, Omega, sigma, w0, gamma, m)
N = numel(t);
h = t(2) - t(1);

Z = zeros(N, 4);
Z(1,:) = z0(:)';

c   = cos(sigma);
w02 = w0^2;

damp = 0;
if gamma ~= 0
    damp = gamma/m;
end

for n = 1:N-1
    z = Z(n,:)';

    k1 = f(z, Omega, c, w02, damp);
    k2 = f(z + 0.5*h*k1, Omega, c, w02, damp);
    k3 = f(z + 0.5*h*k2, Omega, c, w02, damp);
    k4 = f(z + h*k3,     Omega, c, w02, damp);

    Z(n+1,:) = (z + (h/6)*(k1 + 2*k2 + 2*k3 + k4))';
end
end

function dz = f(z, Omega, c, w02, damp)
x  = z(1);  y  = z(2);
vx = z(3);  vy = z(4);

ax =  2*Omega*c*vy - w02*x - damp*vx;
ay = -2*Omega*c*vx - w02*y - damp*vy;

dz = [vx; vy; ax; ay];
end