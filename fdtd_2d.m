%% Define Grid
Nx = 100; Ny = 100;
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
steps = 150;

%% Main Loop
for T = 1 : steps
    % Calculate Dz Field
    for i = 2 : Nx
        for j = 2 : Ny
            Dz(i, j) = Dz(i, j) + 0.5 * (Hy(i, j) - Hy(i-1, j) - Hx(i, j) + Hx(i, j-1));
        end
    end

    % Put a Gaussian pulse in the middle
    pulse = 10 * exp(-0.5*(((t0 - T) / spread)^2));
    Dz(Nx/2, Ny/2) = Dz(Nx/2, Ny/2) + pulse;

    % PECs at the boundaries
    Dz(1, :) = 0;
    Dz(Nx, :) = 0;
    Dz(:, 1) = 0;
    Dz(:, Ny) = 0;

    % Calculate Ez Field
    for i = 1 : Nx
        for j = 1 : Ny
            Ez(i, j) = Ga(i, j) * Dz(i, j);
        end
    end

    % Calculate Hx Field
    for i = 1 : Nx - 1
        for j = 1 : Ny - 1
            Hx(i, j) = Hx(i, j) + 0.5 * (Ez(i, j) - Ez(i, j+1));
        end
    end
    

    % Calculate Hy Field
    for i = 1 : Nx - 1
        for j = 1 : Ny - 1
            Hy(i, j) = Hy(i, j) + 0.5 * (Ez(i+1, j) - Ez(i, j));
        end
    end

    imagesc(xa, ya, Dz');
    colormap 'jet';
    clim([-1 1]);
    colorbar;
    title(['FDTD Simulation (Time step: ', num2str(T), '/', num2str(steps), ')']);
    xlabel('X');
    ylabel('Y');
    drawnow;
end

