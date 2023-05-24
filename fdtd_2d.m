%% Define Grid
Nx = 500; Ny = 500;
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
steps = 500;

% Exporter
videoFile = 'simulation_movie2.mp4';
videoObj = VideoWriter(videoFile, 'MPEG-4');
videoObj.FrameRate = 30;
open(videoObj);

%% Main Loop
for T = 1 : steps
    % Calculate Dz Field
    for i = 2 : Nx
        for j = 2 : Ny
            Dz(i, j) = Dz(i, j) + 0.5 * (Hy(i, j) - Hy(i-1, j) - Hx(i, j) + Hx(i, j-1));
        end
    end

    % Put a Gaussian pulse in the middle
    % pulse = 10 * exp(-0.5*(((t0 - T) / spread)^2));
    % Dz(Nx/2, Ny/2) = Dz(Nx/2, Ny/2) + pulse;

    freq = 1e9;
    lambda = c0/freq;
    radius = lambda;

    center = [Nx/2, Ny/2];
    for i = 1 : Nx
        for j = 1 : Ny
            dist = sqrt(((center(1) - i)*dx)^2 + ((center(2) - j)*dy)^2);
            if abs(dist - radius) <= 3 * dx
                % Calculate angle in degrees (-180 to 180)
                deg = rad2deg(angle((i - center(1)) + 1j * (j - center(2))));
                phase = deg / 2;
                amplitude = -(deg/180)^2 + 1.5;
                sigmoid = 1 ./ (1 + exp(-0.1*(T - 50)));

                % Inject hard source (simulating metal)
                pulse = 10 * sin(2 * pi * freq * T * dt + phase);
                Dz(i, j) = amplitude * sigmoid * pulse;
            end
        end
    end

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

    if mod(T, 1) == 0
        imagesc(xa, ya, Dz');
        colormap('parula');
        clim([-1 1]);
        colorbar;
        title(['FDTD Simulation (Time step: ', num2str(T), '/', num2str(steps), ')']);
        xlabel('X');
        ylabel('Y');
        drawnow;
        writeVideo(videoObj, getframe(gcf));
    end
end

%% Deinitalization
close(videoObj);

