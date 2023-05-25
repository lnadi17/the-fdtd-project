%% Define Grid
Nx = 1000; Ny = 1000;
dx = 0.005;
dy = 0.005;
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

URxx = ones(Nx, Ny);
URyy = ones(Nx, Ny);
mEz1 = ones(Nx, Ny);

% Detector
startstep = 1000;
center = [Nx/2, Ny/2];
detectorvalues = -ones(Nx, Ny);
detectorangles = zeros(Nx, Ny);
for i = 1 : Nx
    for j = 1 : Ny
        sqrdist = (i-center(1))^2 + (j-center(2))^2;
        if abs(sqrdist - 400^2) <= 1000
            detectorvalues(i, j) = 0;
            detectorangles(i, j) = angle((i - center(1)) + 1j * (j - center(2)));
        end
    end
end

% Source
t0 = 20;
spread = 6.0;
beta = 0; % Second source's phase
freq = 1e9;
lambda = c0/freq;
radius = lambda;
center1 = [Nx/2, Ny/2 + 3 * radius/dx];
center2 = [Nx/2, Ny/2 - 3 * radius/dx];

% Simulation
steps = 3000;

% Perfectly Matched Layer
Nx2 = 2*Nx;
Ny2 = 2*Ny;
NPML = [20 21 22 23];

sigx = zeros(Nx2,Ny2);
for nx = 1 : 2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1 : 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end
sigy = zeros(Nx2,Ny2);
for ny = 1 : 2*NPML(3)
    ny1 = 2*NPML(3) - ny + 1;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(3))^3;
end
for ny = 1 : 2*NPML(4)
    ny1 = Ny2 - 2*NPML(4) + ny;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(4))^3;
end

% Update coefficients
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);

mHx0 = (1/dt) + sigHy/(2*e0);
mHx1 = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2 = -c0./URxx./mHx0;
mHx3 = -(c0*dt/e0)*sigHx./URxx./mHx0;

sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);

mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = - c0./URyy./mHy0;
mHy3 = - (c0*dt/e0) * sigHy./URyy./mHy0;

sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);

mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1 ./ mDz0;
mDz2 = c0./mDz0;
mDz4 = - (dt/e0^2)*sigDx.*sigDy./mDz0;

IDz = zeros(Nx, Ny);
ICEy = zeros(Nx, Ny);
ICEx = zeros(Nx, Ny);

CEx = zeros(Nx, Ny);
CEy = zeros(Nx, Ny);
CHz = zeros(Nx, Ny);

% Exporter
videoFile = 'simulation_movie3.mp4';
videoObj = VideoWriter(videoFile, 'MPEG-4');
videoObj.FrameRate = 30;
open(videoObj);

%% Main Loop
for T = 1 : steps
    % Update curls of Ex and Ey
    % Compute CEx
    CEx(1:Nx,1:Ny-1) = (Ez(1:Nx,2:Ny)-Ez(1:Nx,1:Ny-1))/dy;
    CEx(1:Nx,Ny) = -Ez(1:Nx,Ny)/dy;

    % Compute CEy
    CEy(1:Nx-1,1:Ny) = -1*(Ez(2:Nx,1:Ny) - Ez(1:Nx-1,1:Ny))/dx;
    CEy(Nx,1:Ny) = Ez(Nx,1:Ny)/dx;

    % Update H integrations
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;

    % Calculate Hx Field
    Hx = mHx1 .* Hx + mHx2 .* CEx + mHx3 .* ICEx;

    % Calculate Hy Field
    Hy = mHy1 .* Hy + mHy2 .* CEy + mHy3 .* ICEy;

    % Compute curl of H
    CHz(1,1) = (Hy(1,1))/dx - (Hx(1,1))/dy;
    CHz(2:Nx,1) = (Hy(2:Nx,1) - Hy(1:Nx-1,1))/dx - (Hx(2:Nx,1))/dy;
    CHz(1,2:Ny) = (Hy(1,2:Ny))/dx - (Hx(1,2:Ny) - Hx(1,1:Ny-1))/dy;
    CHz(2:Nx,2:Ny) = (Hy(2:Nx,2:Ny) - Hy(1:Nx-1,2:Ny))/dx - ...
        (Hx(2:Nx,2:Ny) - Hx(2:Nx,1:Ny-1))/dy;

    % Update D integrations
    IDz = IDz + Dz;

    % Calculate Dz Field
    Dz = mDz1 .* Dz + mDz2 .* CHz + mDz4 .* IDz;

    % Inject circular source
    for i = 1 : Nx
        for j = 1 : Ny
            dist1 = sqrt(((center1(1) - i)*dx)^2 + ((center1(2) - j)*dy)^2);
            dist2 = sqrt(((center2(1) - i)*dx)^2 + ((center2(2) - j)*dy)^2);
            if abs(dist1 - radius) <= 6 * dx
                % Calculate angle in degrees (-180 to 180)
                deg = rad2deg(angle((i - center1(1)) + 1j * (j - center1(2))));
                phase = deg / 2;
                amplitude = -(deg/180)^2 + 1.5;
                sigmoid = 1 ./ (1 + exp(-0.1*(T - 50)));

                % Inject hard source (simulating metal)
                pulse = 10 * sin(2 * pi * freq * T * dt + phase);
                Dz(i, j) = amplitude * sigmoid * pulse;
            end
            if abs(dist2 - radius) <= 6 * dx
                % Calculate angle in degrees (-180 to 180)
                deg = rad2deg(angle((i - center2(1)) + 1j * (j - center2(2))));
                phase = deg / 2;
                amplitude = -(deg/180)^2 + 1.5;
                sigmoid = 1 ./ (1 + exp(-0.1*(T - 50)));

                % Inject hard source (simulating metal)
                pulse = 10 * sin(2 * pi * freq * T * dt + phase + beta);
                Dz(i, j) = amplitude * sigmoid * pulse;
            end
        end
    end

    % Calculate Ez Field
    Ez = mEz1 .* Dz;

    % Start detector when step count reaches specific value
    if T >= startstep
        for i = 1 : Nx
            for j = 1 : Ny
                if detectorvalues(i, j) >= 0
                    detectorvalues(i, j) = detectorvalues(i, j) + abs(Dz(i, j));
                end
            end
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

%% Post Processing
close(videoObj);

a = [];
detectorvalues = detectorvalues ./ (steps - startstep);
for i = 1 : Nx
    for j = 1 : Ny
        if detectorvalues(i, j) > 0
            a = [a; [detectorangles(i, j), detectorvalues(i, j)]];
        end
    end
end

[G, aa] = findgroups(a(:, 1));
b = [aa, splitapply(@mean, a(:, 2), G)];
polarplot(b(:, 1), b(:, 2))



