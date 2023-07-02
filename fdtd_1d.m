%% Define Grid

% NOTE: Grid resolution should must be sufficient to resolve the shortest
% wavelength. For example, for most metallic surfaces 60-100 cells in one
% wavelength should be enough.

% NOTE: If our simulation has critical dimensions, snapping grid to these
% dimensions should be considered.

Sz = 4; % Physical size along z
Nz = 400; % Number of cells along z
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
dt = dz / 2.0 / c0; % Single time step in seconds
steps = 1000; % Total simulation steps
t = (0:dt:(steps - 1) * dt); % Each time step

% Initialize materials to free space
ER = ones(1, Nz) .* er;
UR = ones(1, Nz) .* ur;
SIGMA = zeros(1, Nz);

% Define geometry
ER(101:200) = ones(1, 100) .* 4;
SIGMA(101:200) = ones(1, 100) .* 0.04;

% Compute update coefficients
EAF = dt * SIGMA ./ (2 * e0 .* ER);
aEy = (1 - EAF) ./ (1 + EAF);
bEy = (c0 * dt) ./ (ER .* (1 + EAF)) ./ dz;

mHx = 0.5 ./ UR;

% Initialize H and E fields to zeros
Hx = zeros(1, Nz);
Ey = zeros(1, Nz);

% Initialize variables for perfect boundary conditions
h1 = 0; h2 = 0;
e1 = 0; e2 = 0;

%% Describe source
fmax = 1e9; % Max frequency we're interested in
tau = 0.5 / fmax; % FWHM (Full Width at Half Maximum)
t0 = 6 * tau; % Pulse offset
th = t + dt + dt/2;

% These sources only work for Ey/Hx mode, and if the source is in the free space
gaussian = exp(-((t - t0) ./ tau ) .^ 2); % Gaussian pulse
gaussianH = -exp(-((th - t0) ./ tau) .^ 2); % Gaussian pulse for directional source H
harmonic = 1 * sin(2 * pi * 700e6 * t); % 700 MHz harmonic source
harmonicH = -1 * sin(2 * pi * 700e6 * th); % 700 MHz harmonic source for directional source H

% Choose the desired pulse for simulation
pulse = harmonic;
pulseH = harmonicH;

% Exporter
videoFile = 'simulation_movie_1d.mp4';
videoObj = VideoWriter(videoFile, 'MPEG-4');
videoObj.FrameRate = 30;
open(videoObj);

%% Main FDTD Loop
jump = 1;
prompt = false;
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
    Ey(1) = aEy(1) * Ey(1) + bEy(1) * (Hx(1) - h2); % Perfect boundary condition
    for nz = 2 : Nz
        Ey(nz) = aEy(nz) * Ey(nz) + bEy(nz) * (Hx(nz) - Hx(nz-1));
    end

    % Inject source
    nzsrc = 5; % Source is in the fifth cell
    
    % Directional source
    Hx(nzsrc - 1) = Hx(nzsrc - 1) - mHx(nzsrc - 1) * pulse(T);
    Ey(nzsrc) = Ey(nzsrc) - bEy(nzsrc) * pulseH(T);

    % Soft source
    % Ey(nzsrc) = Ey(nzsrc) + pulse(T);

    % Hard Source
    % Ey(nzsrc) = pulse(T);

    % Try to put metal in the center
    % Ey(Nz/2) = 0;

    % Visualize E and H (not necessarily every step)
    if mod(T, jump) == 0
        clf;
        hold on;
        plot(Ey, 'LineWidth', 2);
        % plot(Hx, 'LineWidth', 2);
        rectangle('Position', [100, -1.5, 100, 3], 'FaceColor', [0, 0, 1, 0.3], 'EdgeColor', 'none');
        ylim([-1.5, 1.5]);
        xlim([0, Nz]);
        title(['Step ', num2str(T)])
        drawnow;
        writeVideo(videoObj, getframe(gcf));
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

%%

close(videoObj);