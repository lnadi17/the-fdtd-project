%% Define Grid

% NOTE: Grid resolution should must be sufficient to resolve the shortest
% wavelength. For example, for most metallic surfaces 60-100 cells in one 
% wavelength should be enough.

% NOTE: If our simulation has critical dimensions, snapping grid to these
% dimensions should be considered.

Sx = 1; % Physical size along x
Sy = 1; % Physical size along y
Sz = 1; % Physical size along z

Nx = 1000; % Number of cells along x
Ny = 1000; % Number of cells along y
Nz = 1000; % Number of cells along z

dx = Sx/Nx;
dy = Sy/Ny;
dz = Sz/Nz;

%% Define Rectangle Geometry
wx = 0.2;
wy = 0.6;

% Compute position indices
nx = round(wx/dx);
nx1 = 1 + floor((Nx - nx)/2);
nx2 = nx1 + nx - 1;

ny = round(wy/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;

% Create A
A = zeros(Nx,Ny);
A(nx1:nx2,ny1:ny2) = 1;

%% Define Physical Parameters

% Speed of light
c0 = 2.99792458e8;

% Permittivity and permeability of free space
u0 = 4.0 * pi * 1.0e-7;
e0 = 1.0 / (c0 * c0 * u0);
er = 1.0;
ur = 1.0;

% Time step (according to Courant Condition)
dt = min(min(dx, dy), dz) / 2.0 / c0; % single time step in seconds
steps = 500; % total simulation steps
timestamps = (0:dt:(steps - 1) * dt) .* 1e9; % each time step in nanoseconds

% Initialize materials to free space
ER = ones(1, Nz) .* er;
UR = ones(1, Nz) .* ur;

% Compute update coefficients (only for z direction)
mEy = (c0 * dt) ./ ER / dz;
mHx = (c0 * dt) ./ UR / dz;

Hx = zeros(1, Nz);
Ey = zeros(1, Nz);

% For perfect boundary conditions
h1 = 0; h2 = 0;
e1 = 0; e2 = 0;

%% Main FDTD Loop
for T = 1 : steps
    h2 = h1; h1 = Hx(1); % Record boundary Hx that we'll use 2 steps later
    % Update H from E
    for nz = 1 : Nz - 1
        Hx(nz) = Hx(nz) + mHx(nz) * (Ey(nz+1) - Ey(nz));
    end
    Hx(Nz) = Hx(Nz) + mHx(nz) * (e2 - Ey(Nz)); % Perfect boundary condition
    

    % Update E from H
    e2 = e1; e1 = Ey(Nz); % Record boundary Ex that we'll use 2 steps later
    Ey(1) = Ey(1) + mEy(1) * (h2 - Hx(1)); % Perfect boundary condition
    for nz = 2 : Nz
        Ey(nz) = Ey(nz) + mEy(nz) * (Hx(nz) - Hx(nz-1));
    end
end