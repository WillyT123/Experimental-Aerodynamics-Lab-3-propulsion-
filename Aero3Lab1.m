clear all
clc
tic
instr_data = cell(2,10);
volts = cell(10,3);
% Intstrumentation Data Layout
% P_transducer(psf) |   P_pitot(Pa) |   T_acqired(C)    |   T_wire(C)   |
% T_ref(C)  |   E_acqired(V DC)
rho = 1.14; %kg/m^3
% E_corrected = sqrt((T_wire - T_ref)/(T_wire-T_acq))*E_acq
for h = 1:10
    Negfiles = 'instrumentation_data_Q_NEG%d.dat';
    Posfiles = 'instrumentation_data_Q_POS%d.dat';
instr_data{1,h} = load(sprintf(Negfiles,h));
instr_data{2,h} = load(sprintf(Posfiles,h));
end

E_n = zeros(3,10);
E_p = zeros(3,10);
for h = 1:10
E_n(:,h) = sqrt((instr_data{1,h}(:,5)-instr_data{1,h}(:,6))/(instr_data{1,h}(:,5)-instr_data{1,h}(:,4)))*instr_data{1,h}(:,7);
E_p(:,h) = sqrt((instr_data{2,h}(:,5)-instr_data{2,h}(:,6))/(instr_data{2,h}(:,5)-instr_data{2,h}(:,4)))*instr_data{2,h}(:,7);
end

V_pit = zeros(3,20);
for h = 1:10
V_pit(:,h) = abs((2*(instr_data{1,h}(:,3))/rho).^(1/2));
V_pit(:,h+10) = abs(sqrt(2*(instr_data{2,h}(:,3))/rho));
end



%%
figure
hold on

for h = 1:3
plot(E_n(h,:),V_pit(h,1:10),'*')
plot(E_p(h,:),V_pit(h,11:20),'o')
end
ylabel('Pitot Velocity (m/s)')
xlabel('CTA Voltage (V DC)')
title('Voltage vs Velocity')
legend('Neg curve','Pos curve')
hold off

degree = 2;

phit = zeros(6,degree+1);
f1  = zeros(6,10);
figure
hold on

for h = 1:3
phit(h,:) = polyfit(E_n(h,:),V_pit(h,1:10),degree)
phit(h+3,:) = polyfit(E_p(h,:),V_pit(h,11:20),degree)
f1(h,:) = polyval(phit(h,:),E_n(h,:));
f1(h+3,:) = polyval(phit(h+3,:),E_p(h,:));
plot(E_n(h,:),f1(h,:))
plot(E_p(h,:),f1(h+3,:))
end
fit = mean(phit,1);
ylabel('Pitot Velocity (m/s)')
xlabel('CTA Voltage (V DC)')
title('Fitted Voltage Curves')
legend('Neg curve','Pos curve')
hold off

%%
mu = 1.962e-5;
d = 0.1016;

tu = zeros(10,1);
figure
v_avg = zeros(10,1);
for h = 1:10
    %for y = 1:3
   Filename = 'voltage_time_history_Q_POS%d_position_%d.dat';     
        volts{h,2} = load(sprintf(Filename,h,2));
    %end
    vel = f1(5,1)*(volts{h,2}(:,2)).^3 +f1(5,2)*(volts{h,2}(:,2)).^2 + f1(5,3)*(volts{h,2}(:,2))+f1(5,4);
    v_avg(h) = mean(vel);
    tu(h) = (std(vel)/mean(vel))*100; 
end

plot(v_avg,tu)
xlabel('Average Velocity (m/s)')
ylabel('Turbulence %')
title('Turbulence Curve')
hold off
%%

cd lab2_cta_vortex

mat = dir('*.dat');

p = single(0.5);
lab2 = cell(17,2);
for h = 1:length(mat)
 files  = mat(h).name;
if files(28)=='l' || files(29)=='l'
    loading = importdata(files);         %Laminar sorted into col 1 of cell array, turbulent in col 2
    lab2{h,1} = loading(:,2);            % This convention is maintained for the remainder
elseif files(28)=='t' || files(29)=='t'
    loading = importdata(files);
    lab2{h,2} = loading(:,2);
end
end
cd C:\Users\Travis\Documents\MATLAB\PlaneStuff\Aero3
%%

r = 1;
d = 0.1016;

L = 98304;
Fs = 10000;
f = Fs*(0:(L/2))/L;
i = find(f>3);
j = find(f>500);

Hz = zeros(17,2);
for t = 1:17
Fourier(:,1) = fft(lab2{r,1});      %
Fourier(:,2) = fft(lab2{r+1,2});
P2 = abs(Fourier/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);


[a,l] = max(P1(i(1):j(1),1));      %Laminar indicies stored in l
[b,u] = max(P1(i(1):j(1),2));      %Turbulent indicies stored in u

Hz(t,1) = f(l);
Hz(t,2) = f(u);

% plot(f(i(1):j(1)),P1(i(1):j(1))) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

r = r + 2;
end

press = [0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 6 7 8 9 10 11 12]'*47.88; %Hardcoded dyn pressures

V_trans = (sqrt(2*press/rho));    %Resulting hardcoded velocity and Re
Re = rho*V_trans*d/mu;

strouhal(1:17,1) = Hz(:,1)*d./V_trans;     %Laminar Strouhal values
strouhal(1:17,2) = Hz(:,2)*d./V_trans;     %Turbulent Strouhal values
figure
hold on
plot(Re,strouhal(:,1))
plot(Re,strouhal(:,2))
title('Strouhal vs Re analysis')
xlabel('Reynolds #')
ylabel('Strouhal #')
legend('Laminar data','Turbulent data')

%%
% Create BS points within the figure's limits
% (don't try to plot something at x=2^30 when the xlim is at x=5...)
% [X, Y] = pol2cart([0:0.01:(2.*pi)], cos([0:0.01:(2.*pi)])+cos(2.*[0:0.01:(2.*pi)]));
% X = X+abs(min(X));
% Y = Y+abs(min(Y));
% X = X.*(10^6-10^5)./(max(X)-min(X));
% Y = Y.*0.5./(max(Y)-min(Y));
% Less idiotic data for you fools out there:
% X = [0, linspace(10^5, 8*10^5, 20), 10^6];
% Y = [0, 0.1+0.3.*rand(1, 20), 0.5];


X = Re;
Y = strouhal(:,2);

 %MLG 3000 patent pending code:
% Use the code from here, with your own X and Y data, saved as 'X' and 'Y'
p_X = [107, 890]; % Image's X limits, in pixels, of the plot's position
p_Y = [88, 746]; % Same for Y
X_lim = [4*10^4, 4*10^5]; % X axis limits on the figure you want to plot
Y_lim = [0, 0.5];

X_plot = p_X(1) + X.*diff(p_X)./diff(X_lim); % Magic to scale things onto plot
Y_plot = p_Y(1) + Y.*diff(p_Y)./diff(Y_lim);

figure('Name', 'MLG 3000 patent pending code', 'NumberTitle', 'off');
Fig = flip(imread('MLG3000_plot.JPG'), 1);
[stdby1, stdby2, ~] = size(Fig);
imshow(Fig);
set(gca, 'YDir', 'normal')
line(X_plot, Y_plot, 'color', 'r', 'lineWidth', 2);
print(gcf, 'MLG3000_output.png', '-dpng', '-r300');

toc



