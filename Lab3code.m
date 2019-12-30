clc
clear all
close all
tic
data = load('Lab-3_prop_data_all_sessions.dat');


% No. of Props | Prop dia (in) | Prop Pitch (in) | Throttle | Q (psf) | V (Volts) | I (Amps) 
%                              | Thrust (lb) | Torque (in-lb) | RPM | Material (1 - APC; 2- Wood)

rho = 1.18; %kg/m^3

%co-plot the following curves once for each prop
%Ct v J, Cq v J, Cp v J, eta v J (7 total plots)

%co-Plot varying pitch, diameter, and blade count v J once each for
%variations in diam., pitch, and blade count
data(:,2) = data(:,2)*0.0254; %in to m
data(:,3) = data(:,3)*0.0254; %in to m
data(:,5) = data(:,5)*47.77; %psf to pa
data(:,8) = data(:,8)*4.4482; %lbs to newtons
data(:,9) = data(:,9)*0.113; %in-lbs to N-m
data(:,10) = data(:,10)./60; % rot/m to rot/s
data(199,:)= zeros(1,11);   %add line of zeros to data to trigger if-statement last time
%data(data<0) = NaN;


Vinf = zeros(32,7);
J = zeros(32,7);
Ct = zeros(32,7);
Cq = zeros(32,7);           %pre-allocating
P = zeros(32,7);
Cp = zeros(32,7);
eta = zeros(32,7);

r = 1;
c = 1;                      %initializing matrix references
ref1 = 1;

d = data(1,2);
p = data(1,3);

for i = 1:199
   if d~=data(i,2) || p~=data(i,3) 
       
       ref2 = i-1;
       
       Vinf(1:r-1,c) = sqrt(2*data(ref1:ref2,5)/rho);
       J(1:r-1,c) = Vinf(1:r-1,c)./(data(ref1:ref2,10).*data(ref1:ref2,2));
       Ct(1:r-1,c) = data(ref1:ref2,8)./(rho*data(ref1:ref2,10).^2.*data(ref1:ref2,2).^4);
       Cq(1:r-1,c) = data(ref1:ref2,9)./(rho*data(ref1:ref2,10).^2.*data(ref1:ref2,2).^5);
       P(1:r-1,c) = data(ref1:ref2,6).*data(ref1:ref2,7);
       Cp(1:r-1,c) = P(1:r-1,c)./(rho*data(ref1:ref2,10).^3.*data(ref1:ref2,2).^5);      
       eta(1:r-1,c) = Ct(1:r-1,c).*J(1:r-1,c)./(2*pi.*Cq(1:r-1,c));
        
       r = 2;
       c = c+1;
       ref1 = i;
   else
        r = r+1;
   end
        p = data(i,3);
        d = data(i,2);   
end

%  Propellers are like wings but they spin. 
% Spinnny spinny wing thing
% Don’t touch spinny wing 
% You lose touchy thing

Cq = 10*Cq;
Cq(Cq==0) = nan;
Ct(Ct==0) = nan;
Cp(Cp==0) = nan;
eta(eta>1) =nan;
eta(eta<0) = nan;

di = single([16 15 10 10 13 10 10.5]);
pi = [8 8 7 6 8 8 8];

for j = 1:7

avr = J(:,j);
thr = Ct(:,j);
tor = Cq(:,j);

airbrake = min(avr(thr<0));
windmill = min(avr(tor<0 & thr<0));


figure(j)
hold on
grid on
plot(J(:,j),Ct(:,j),'*')
plot(J(:,j),Cq(:,j),'+')
plot(J(:,j),Cp(:,j),'s')
plot(J(:,j),eta(:,j),'o')
xline(airbrake,'r-.');
xline(windmill,'--');
title(sprintf('%2.1fx%2.1f Propeller',di(j),pi(j)))
xlabel('Advanced Ratio')
legend('Ct v J','Cq v J','Cp v J','eta v J','Airbrake Regime','Windmill Regime')
hold off
end
toc
%%
%Diameters (props 1,2,5,6)

figure
grid on
hold on
for s = [1 2 5 6]
plot(J(:,s),Ct(:,s),'*')
end
title('Ct v J, Diam.')
xlabel('Advanced Ratio')
ylabel('Coeff. of thrust')
legend('16x8','15x8','13x8','10x8')
hold off

figure
grid on
hold on
for s = [1 2 5 6]
plot(J(:,s),Cq(:,s),'+')
end
title('Cq v J, Diam.')
xlabel('Advanced Ratio')
ylabel('Coeff. of torque')
legend('16x8','15x8','13x8','10x8')
hold off


figure
grid on
hold on
for s = [1 2 5 6]
plot(J(:,s),Cp(:,s),'s')
end
title('Cp v J, Diam.')
xlabel('Advanced Ratio')
ylabel('Coeff. of power')
legend('16x8','15x8','13x8','10x8')
hold off


figure
hold on
grid on
for s = [1 2 5 6]
plot(J(:,s),eta(:,s),'o')
end
title('eta v J, Diam.')
xlabel('Advanced Ratio')
ylabel('efficiency')
legend('16x8','15x8','13x8','10x8')
hold off

%%
%Pitch

figure
hold on
grid on
for s = [3 4 6]
plot(J(:,s),Ct(:,s),'*')
end
title('Ct v J, Pitch')
xlabel('Advanced Ratio')
ylabel('Coeff. of thrust')
legend('10x7','10x6','10x8')
hold off

figure
hold on
grid on
for s = [3 4 6]
plot(J(:,s),Cp(:,s),'s')
end
title('Cp v J, Pitch')
xlabel('Advanced Ratio')
ylabel('Coeff. of power')
legend('10x7','10x6','10x8')
hold off

figure
grid on
hold on
for s = [3 4 6]
plot(J(:,s),Cq(:,s),'+')
end
title('Cq v J, Pitch')
xlabel('Advanced Ratio')
ylabel('Coeff. of torque')
legend('10x7','10x6','10x8')
hold off

figure
grid on
hold on
for s = [3 4 6]
plot(J(:,s),eta(:,s),'o')
end
title('eta v J, Pitch')
xlabel('Advanced Ratio')
ylabel('efficiency')
legend('10x7','10x6','10x8')
hold off

%%

%blade number

figure
grid on
hold on
for s = [6 7]
plot(J(:,s),eta(:,s),'o')
end
title('eta v J, Blade #')
xlabel('Advanced Ratio')
ylabel('efficiency')
legend('2 blades','4 blades')
hold off


figure
grid on
hold on
for s = [6 7]
plot(J(:,s),Cp(:,s),'s')
end
title('Cp v J, Blade #')
xlabel('Advanced Ratio')
ylabel('Coeff. of power')
legend('2 blades','4 blades')
hold off


figure
grid on
hold on
for s = [6 7]
plot(J(:,s),Ct(:,s),'*')
end
title('Ct v J, Blade #')
xlabel('Advanced Ratio')
ylabel('Coeff. of thrust')
legend('2 blades','4 blades')
hold off

figure
grid on
hold on
for s = [6 7]
plot(J(:,s),Cq(:,s),'+')
end
title('Cq v J, Blade #')
xlabel('Advanced Ratio')
ylabel('Coeff. of torque')
legend('2 blades','4 blades')
hold off