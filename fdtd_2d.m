%% Define Grid
Nx = 60; Ny = 60;
dx = 0.01;
dy = 0.01;
Sx = Nx*dx;
Sy = Ny*dy;

xa = (1:Nx)*dx;
ya = (1:Ny)*dy;
[Y,X] = meshgrid(ya,xa);

%% Constants
c0 = 2.99792458e8;
dt = dx/2/c0;
u0 = 4.0 * pi * 1.0e-7;
e0 = 1.0 / (c0 * c0 * u0);

%% Initialization
% Fields
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
Ez = zeros(Nx, Ny);
Dz = zeros(Nx, Ny);
Ga = ones(Nx, Ny);

% Source
t0 = 20;
spread = 6.0;

% Simulation
steps = 200;

%% Main Loop
for T = 1 : steps
    % Calculate Dz Field
    for i = 2 : Nx - 1
        for j = 2 : Ny - 1
            Dz(i, j) = Dz(i, j) + 0.5 * (Hy(i, j) - Hy(i-1, j) - Hx(i, j) + Hx(i, j-1));
        end
        Dz(i, 1) = Dz(i, 1) + 0.5 * (Hy(i, 1) - Hy(i-1, 1) - Hx(i, 1) + 0);
    end
    Dz(1, 1) = Dz(1, 1) + 0.5 * (Hy(1, 1) - 0 - Hx(1, 1) + 0);

    % Put a Gaussian pulse in the middle
    pulse = exp(-0.5*(((t0 - T) / spread)^2));
    Dz(Nx/2, Ny/2) = Dz(Nx/2, Ny/2) + pulse;

    % Calculate Ez Field
    for i = 1 : Nx - 1
        for j = 1 : Ny - 1
            Ez(i, j) = Ga(i, j) * Dz(i, j);
        end
    end

    % Calculate Hx Field
    for i = 1 : Nx - 1
        for j = 1 : Ny - 2
            Hx(i, j) = Hx(i, j) + 0.5 * (Ez(i, j) - Ez(i, j+1));
        end
        Hx(i, Ny-1) = Hx(i, Ny-1) + 0.5 * (Ez(i, Ny-1) - 0);
    end
    

    % Calculate Hy Field
    for i = 1 : Nx - 2
        for j = 1 : Ny - 1
            Hy(i, j) = Hy(i, j) + 0.5 * (Ez(i+1, j) - Ez(i, j));
        end
    end
    Hy(Nx-1, Ny-1) = Hy(Nx-1, Ny-1) + 0.5 * (0 - Ez(Nx-1, Ny-1));

    imagesc(xa, ya, Dz');
    pause(0.1);
end

