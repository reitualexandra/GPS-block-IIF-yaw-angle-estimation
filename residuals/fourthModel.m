clear all
close all
clc

Ne = 300; % nr of epochs (should be lower than 361 for 2 hours maneuvers and 30s sampled data!)
Ns = 5; % nr of stations
stationsUsed = [24, 48, 59, 64, 78]; % indices of all stations which are measuring at every epoch - has to be selected manually...

residualFile = "YAW210590M1G27.res";
matData = dlmread(residualFile);

prnSat = matData(1, 2);
stations = matData(:, 1);
epochs = matData(:, 3);
residuals = matData(:, 4);
azimuthAngle = matData(:, 5);
nadirAngle = matData(:, 6) + 90;

lightSpeed = 299792458;

matDataReduced = zeros(Ns, Ne);
epochsReduced = zeros(1, Ne);
azimuthReduced = zeros(Ns, Ne);
for i=1:Ns % length(stationsUsed)
    res = residuals(stations==stationsUsed(i));
    matDataReduced(i,:) = res(1:Ne);
    
    epo = epochs(stations==stationsUsed(i));
    epochsReduced(i,:) = epo(1:Ne);
    
    azi = azimuthAngle(stations==stationsUsed(i));
    azimuthReduced(i,:) = azi(1:Ne);
    
    nad = nadirAngle(stations==stationsUsed(i));
    nadirReduced(i,:) = nad(1:Ne);
end

residualsFinal = [];
nadirFinal = [];
azimuthFinal = [];
for i=1:Ne
    residualsFinal = [residualsFinal; matDataReduced(:, i)];
    nadirFinal = [nadirFinal; nadirReduced(:, i)];
    azimuthFinal = [azimuthFinal; azimuthReduced(:, i)];
end

%% Plot residuals

figure()
scatter(epochsReduced, matDataReduced, '.')
datetick('x','HH:MM')
grid on
title(strcat("Residuals for DOY 059, PRN ", num2str(prnSat)))
xlabel("Time [HH:MM]")
ylabel("Residuals [m]")

%% Create design matrix for Ne epochs and Ns stations
A = [];
B = [];
C = [];
D = [];
for i=1:Ne
    j = Ns*(i-1) + 1;
    nk = deg2rad(nadirFinal(j:j+Ns-1));
    ak = deg2rad(azimuthFinal(j:j+Ns-1));
    
    a = [ sin(nk).*sin(ak) sin(nk).*cos(ak) ];
    A = blkdiag(A, a);
    
end
%%
for i=1:Ne
    x = (A'*A)\A'*residualsFinal;
    k = 2*(i-1) + 1;
    b = [2*x(k) 2*x(k+1)];
    B = blkdiag(B, b);
    
    c = zeros(Ns, Ne);
    c(:, i) = ones(Ns, 1);
    C = [C; c];
    
    d = eye(Ns);
    D = [D; d];
end

%% Create grand design matrix and get solution vector

% A = [A C D; B zeros(Ne, Ne+Ns)]; % with ambiguities estimation - for
% phase residuals
A = [A C ; B zeros(Ne, Ne)]; % without ambiguities estimation - for code residuals
x = (A'*A)\A'*[residualsFinal; 2*(0.394^2)*ones(Ne,1)];

%% Compute yaw angle at each epoch

% ambiguities = x(3*Ne+1:end); % commend when dealing with code estimations
clockCorr = x(2*Ne+1:3*Ne)/lightSpeed;
yaw = zeros(Ne, 1);
for i=1:Ne
    w = x(2*i-1:2*i);
    yaw(i) = atan2(w(2), w(1));
    
end


%% Plot yaw

figure()
scatter(epochsReduced, rad2deg(yaw), '.')
datetick('x','HH:MM')
grid on
title(strcat("Yaw for DOY 059, PRN ", num2str(prnSat)))
xlabel("Time [HH:MM]")
ylabel("Yaw [°]")

