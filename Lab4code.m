clear
clc
clf
close all

% Column definitions in Lab-3_*-Session-Engine-Data.dat: RPM 1 | RPM 2 | RPM 3 | RPM 4 | RPM 5 | Force 1 (gms) |
% Force 2 (gms) | Force 3 (gms) | Force 4 (gms) |Force 5 (gms) | Force Correction (gms) | Time for 10cc fuel consumption (seconds)
% 
% Notes:
% - The RPM measurement is for two blades. You will have to halve the value to get the actual RPM measurement.
% - Average out the five RPM and force measurements.
% - Your torque and power curve trends must be similar to the trends presented during the theory session.
% - Your final data should be in S.I. units.
% - Refer to the lab handout for all constants required for your calculations.

%
A = importdata('C:\Users\Travis\Documents\MATLAB\PlaneStuff\Aero3\Aero3_lab4\Lab4\Lab4.dat'); %imprt the data
%A = importdata('/MATLAB Drive/MAE 451/lab_4.dat');

avg_rpm = (mean(A(:,1:5),2)/2)'; %account for two blades
avg_force = (mean(A(:,6:10)-A(:,11),2).*0.0098)'; %convert the force to newtons from grams

%% Define Variables
avg_rps = avg_rpm/60; %conversion to s
r = .445; %in meters
rho = 862.75; % fuel density in kg/m^3
Qhv = 11300; %heating value of fuel in kJ/kg
nc = .99; %99 percent combustion efficiency

%% calculate stuff
tor = avg_force*r;
pow = 2.*pi.*avg_rps.*tor/1000;
mdot = rho* 1e-5./A(:,12);
eta = pow./(mdot'.*Qhv.*nc);


[a,b] = max(pow);
Epeak = avg_rps(b)
%% Plotting
figure(1)
hold on
scatter(avg_rpm, tor);
xlabel('RPM');
ylabel('Torque (N-M)');
coeff = polyfit(avg_rpm, tor,2);
new_x = linspace(min(avg_rpm), max(avg_rpm), 1000);
new_y = polyval(coeff, new_x);
plot(new_x,new_y);
hold off

figure(2)
hold on 
scatter(avg_rpm, pow);
xlabel('RPM');
ylabel('Engine Power kW');
coeff2 = (polyfit(avg_rpm, pow ,2));
new_x2 = linspace(min(avg_rpm), max(avg_rpm), 1000);
new_y2 = polyval(coeff2, new_x2);
plot(new_x2,new_y2);
hold off

figure(3)
hold on
scatter(avg_rpm, eta);
xlabel('RPM');
ylabel('Engine Efficiency');
coeff3 = polyfit(avg_rpm, eta ,2);
new_x3 = linspace(min(avg_rpm), max(avg_rpm), 1000);
new_y3 = polyval(coeff3, new_x3);
plot(new_x3,new_y3);
hold off


