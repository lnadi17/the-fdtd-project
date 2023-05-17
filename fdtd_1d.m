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

%% Define Physical Parameters

% Speed of light
c0 = 2.99792458e8;

% Permittivity and permeability of free space
u0 = 4.0 * pi * 1.0e-7;
e0 = 1.0 / (c0 * c0 * u0);
er = 1.0;
ur = 1.0;

% Time step (according to Courant Condition)
dt = min(min(dx, dy), dz) / 2.0 / c0; % Single time step in seconds
steps = 3000; % Total simulation steps
t = (0:dt:(steps - 1) * dt); % Each time step

% Initialize materials to free space
ER = ones(1, Nz) .* er;
UR = ones(1, Nz) .* ur;

% Compute update coefficients (only for z direction)
mEy = (c0 * dt) ./ ER / dz;
mHx = (c0 * dt) ./ UR / dz;

% Initialize H and E fields to zeros
Hx = zeros(1, Nz);
Ey = zeros(1, Nz);

% Initialize variables for perfect boundary conditions
h1 = 0; h2 = 0;
e1 = 0; e2 = 0;

%% Describe source
fmax = 1e10; % Max frequency we're interested in
tau = 0.5 / fmax; % tau, aka FWHM (Full Width at Half Maximum)
t0 = 6 * tau; % Pulse offset

gaussian = exp(-((t - t0) ./ tau ) .^ 2); % Gaussian pulse
harmonic = 0.5 * sin(2 * pi * 1e9 * t); % 10 GHz harmonic source

% Choose the desired pulse for simulation
pulse = gaussian;

%% Main FDTD Loop
jump = 20;
prompt = true;
figure;

for T = 1 : steps
    h2 = h1; h1 = Hx(1); % Record boundary H that we'll use 2 steps later
    % Update H from E
    for nz = 1 : Nz - 1
        Hx(nz) = Hx(nz) + mHx(nz) * (Ey(nz+1) - Ey(nz));
    end
    Hx(Nz) = Hx(Nz) + mHx(nz) * (e2 - Ey(Nz)); % Perfect boundary condition

    % Update E from H
    e2 = e1; e1 = Ey(Nz); % Record boundary E that we'll use 2 steps later
    Ey(1) = Ey(1) + mEy(1) * (Hx(1) - h2); % Perfect boundary condition
    for nz = 2 : Nz
        Ey(nz) = Ey(nz) + mEy(nz) * (Hx(nz) - Hx(nz-1));
    end

    % Inject source
    nzsrc = Nz / 4; % Source is at the quarter of the length
    
    % Soft source
    Ey(nzsrc) = Ey(nzsrc) + pulse(T);

    % Hard Source
    % Ey(nzsrc) = pulse(T);

    % Visualize E and H (not necessarily every step)
    if mod(T, jump) == 0
        clf;
        hold on;
        plot(Ey);
        plot(Hx);
        ylim([-1.5, 1.5]);
        xlim([0, Nz]);
        title(['Step ', num2str(T)])
        pause(0.05);
    end

    if T == 500 && prompt == true
        % Prompt the user to continue or stop
        answer = questdlg('Do you want to continue?', 'Continue Simulation', 'Continue', 'Stop', 'Continue');

        if strcmp(answer, 'Stop')
            disp('Simulation stopped by the user.');
            break; % Exit the loop
        end
    end
end