close all
clc
clear
tic
% Column definitions in motor*-run*.dat: RPM | Torque (N-mm) | Voltage (V) | Current (I)
% 
% Notes:
% - You will be presenting a total of 4 plots showing the results for all six data sets.
%      torque, Power in, Power out, efficiency vs rpm

%P_in = I*V
%P_out = torque*rpm
%eff = P_out/P_in

% - Plot scatter plots instead of line plots as your data might have outliers.
% - Your final data should be in S.I. units.


%
FileName = 'motor%d-run%d.dat';

A = cell(2,3);
P_out = cell(2,3);
P_in = cell(2,3);
eff = cell(2,3);

figure
hold on

for a = 1:2
    for b = 1:3
A{a,b} = load(sprintf(FileName,a,b)); %load the data
rads = A{a,b}(:,1)*pi/30;
A{a,b}(:,2) = A{a,b}(:,2)*0.001;

P_out{a,b} = rads.*A{a,b}(:,2);
P_in{a,b} = A{a,b}(:,3).*A{a,b}(:,4);

eff{a,b} = 100*P_out{a,b}(:,1)./P_in{a,b}(:,1);

figName = 'Motor %d, run &d';


subplot(2,2,1)
plot(A{a,b}(:,1),A{a,b}(:,2))
ylabel('Torque (N-m)')
xlabel('RPMs')
legend('Motor 1 run 1','Motor 1 run 2','Motor 1 run 3','Motor 2 run 1','Motor 2 run 2','Motor 2 run 3')
grid on
hold on
subplot(2,2,2)
plot(A{a,b}(:,1),P_in{a,b}(:,1))
ylabel('Power In (w)')
xlabel('RPMs')
legend('Motor 1 run 1','Motor 1 run 2','Motor 1 run 3','Motor 2 run 1','Motor 2 run 2','Motor 2 run 3')
grid on
hold on
subplot(2,2,3)
plot(A{a,b}(:,1),P_out{a,b}(:,1))
ylabel('Power Out (w)')
xlabel('RPMs')
legend('Motor 1 run 1','Motor 1 run 2','Motor 1 run 3','Motor 2 run 1','Motor 2 run 2','Motor 2 run 3')
grid on
hold on
subplot(2,2,4)
plot(A{a,b}(:,1),eff{a,b}(:,1))
ylabel('% Efficiency')
xlabel('RPMs')
legend('Motor 1 run 1','Motor 1 run 2','Motor 1 run 3','Motor 2 run 1','Motor 2 run 2','Motor 2 run 3')
grid on
hold on
    end
end

%%

[M1,I1] = max(P_out{1,3});
[M2,I2] = max(P_out{2,3});

PP_rpm = [A{1,3}(I1,1),A{2,3}(I2,1)]/(2*pi)

toc


