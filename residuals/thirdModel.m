clear all
close all
clc


residualFile = "YAW210590M1G27.res";
matData = dlmread(residualFile);

prnSat = matData(1, 2);
epochs = matData(:, 3);
residuals = matData(:, 4);
azimuthAngle = matData(:, 5);
nadirAngle = matData(:, 6) + 90;

c = 299792458;

%% Plot residuals

figure()
scatter(epochs, residuals, '.')
datetick('x','HH:MM')
grid on
title(strcat("Residuals for DOY 059, PRN ", num2str(prnSat)))
xlabel("Time [HH:MM]")
ylabel("Residuals [m]")

%% Solve for yaw angle

uniqueEpochs = sort(unique(epochs));
yaw = zeros(length(uniqueEpochs), 1);
ts = zeros(length(uniqueEpochs), 1);

for i=1:length(uniqueEpochs)
    dataIndices = find(epochs == uniqueEpochs(i));
    b = residuals(dataIndices);
    A = zeros(length(dataIndices), 2);
    nk = nadirAngle(dataIndices);
    ak = azimuthAngle(dataIndices);
    A(:,1) = sin(deg2rad(nk)).*sin(deg2rad(ak));
    A(:,2) = sin(deg2rad(nk)).*cos(deg2rad(ak));
    
    x = (A'*A)\A'*b;
    yaw(i) = atan2(x(2), x(1));
    
    A = [A; 2*0.394*cos(yaw(i)) 2*0.394*sin(yaw(i))];
    A = [A [ones(length(dataIndices) , 1); 0]];
    b = [b; 2*(0.394^2)];
    x = (A'*A)\A'*b;
    yaw(i) = atan2(x(2), x(1));
    ts(i) = x(3)/c;
end

%% Plot yaw angle

figure()
scatter(uniqueEpochs, rad2deg(yaw), '.')
datetick('x','HH:MM')
grid on
title(strcat("Estimated yaw angle for DOY 059, PRN ", num2str(prnSat)))
xlabel("Time [HH:MM]")
ylabel("Yaw [°]")

%% Plot satellite clock correction

figure()
scatter(uniqueEpochs, ts, '.')
datetick('x','HH:MM')
grid on
title(strcat("Estimated satellite clock correction for DOY 059, PRN ", num2str(prnSat)))
xlabel("Time [HH:MM]")
ylabel("dt [s]")