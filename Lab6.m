% No. of Props | Prop dia (in) | Prop Pitch (in) | Throttle | Q (psf) |V (Volts) | I (Amps) | Thrust (lb) | Torque (in-lb) | RPM | Material (1 - APC; 2- Wood)
close all
clear all
clc
% Notes:
% - Assume density = 1.18 kg/m^3 for your calculations. 
% - If trust/torque is negative, set the efficiency value to 0. 
% 
A = importdata('C:\Users\Travis\Documents\MATLAB\PlaneStuff\Aero3\Aero3_lab4\Lab4\Lab4.dat'); %imprt the data
%A = importdata('/MATLAB Drive/MAE 451/lab_4.dat');

avg_rpm = (mean(A(:,1:5),2)/2)'; %account for two blades
avg_force = (mean(A(:,6:10)-A(:,11),2).*0.0098)'; %convert the force to newtons from grams
avg_rps = avg_rpm/60; %conversion to s
r = .445; %in meters
rho = 862.75; % fuel density in kg/m^3
Qhv = 11300; %heating value of fuel in kJ/kg
nc = .99; %99 percent combustion efficiency
tor = avg_force*r;
pow = 2.*pi.*avg_rps.*tor/1000;
mdot = rho* 1e-5./A(:,12);
eta_thermal = pow./(mdot'.*Qhv.*nc);

[a,b] = max(pow);
Epeak = avg_rps(b);

%%
cd Propeller Data-Cleaned-20191123
data = load('Lab-3_prop_data_all_sessions_clean.dat');
cd C:\Users\Travis\Documents\MATLAB\PlaneStuff\Aero3\Aero3_lab6
rho = 1.18; %kg/m^3`
g = 9.8;
Qhv = 11.3e6;
eta_c = 0.99;
%co-plot the following curves once for each prop
%Ct v J, Cq v J, Cp v J, eta v J (7 total plots)

%co-Plot varying pitch, diameter, and blade count v J once each for
%variations in diam., pitch, and blade count
data(:,2:3) = data(:,2:3)*0.0254; %in to m
data(:,5) = data(:,5)*47.77; %psf to pa
data(:,8) = data(:,8)*4.4482; %lbs to newtons
data(:,9) = data(:,9)*0.113; %in-lbs to N-m
data(:,10) = data(:,10)./60; % rot/m to rot/s
data(165,:)= zeros(1,11);   %add line of zeros to data to trigger if-statement last time
%data(data<0) = NaN;


Vinf = zeros(29,7);
J = zeros(29,7);
Ct = zeros(29,7);
Cq = zeros(29,7);           %pre-allocating
P_in = zeros(29,7);
Cp = zeros(29,7);
eta_curve = zeros(29,7);
Peng = zeros(1,7);

r = 1;
c = 1;                      %initializing matrix references
ref1 = 1;

d = data(1,2);
p = data(1,3);

for i = 1:165
   if d~=data(i,2) || p~=data(i,3) 
       
       ref2 = i-1;
       
       Vinf(1:r-1,c) = sqrt(2*data(ref1:ref2,5)/rho);
       J(1:r-1,c) = Vinf(1:r-1,c)./(data(ref1:ref2,10).*data(ref1:ref2,2));
       Ct(1:r-1,c) = data(ref1:ref2,8)./(rho*data(ref1:ref2,10).^2.*data(ref1:ref2,2).^4);
       Cq(1:r-1,c) = data(ref1:ref2,9)./(rho*data(ref1:ref2,10).^2.*data(ref1:ref2,2).^5);
       P_in(1:r-1,c) = data(ref1:ref2,6).*data(ref1:ref2,7);
       Cp(1:r-1,c) = P_in(1:r-1,c)./(rho*data(ref1:ref2,10).^3.*data(ref1:ref2,2).^5);      
       eta_prop(1:r-1,c) = Ct(1:r-1,c).*J(1:r-1,c)./(2*pi.*Cq(1:r-1,c));
      
       tork(1:r-1,c) = data(ref1:ref2,9);
       spinny(1:r-1,c) = data(ref1:ref2,10);
       Diam(c,1) = data(ref2,2);
       
       r = 2;
       c = c+1;
       ref1 = i;
   else
       r=r+1;
   end
        p = data(i,3);
        d = data(i,2);   
end


%%

%From lab 5 data

peaks = [932.1, 689.2 Epeak]; %rot/s

%Added Code:
% [M1,I1] = max(P_out{1,3});
% [M2,I2] = max(P_out{2,3});
% 
% PP_rpm = [A{1,3}(I1,1),A{2,3}(I2,1)]/(2*pi)


eff_fit = zeros(7,4);
% for a = 1:7
%    eff_fit(a,:) = polyfit(J(:,a),eta(:,a),3);
%   
% end
g = 1;
% for M=1:3
% for af = 1:2

for af = 1:2
    for M=1:3
if af == 1
    OV = linspace(0,40,29);
    Cd0 = 0.01;
    S = 0.5;
    W1 = 3;
    e = 0.8;
    AR = 3;
    t_f = 15;
else
    OV = linspace(0,60,29);
    Cd0 = 0.018;
    S = 0.8;
    W1 = 10;
    e = 0.6;
    AR =5;
    t_f = 30;
end

for c = 1:7

mdot = P_in(:,c)./(Qhv*eta_c);

W2 = W1*g - (t_f*mdot*g)';

CL = 2*W1./(rho*OV.^2*S);

CD = Cd0 + CL.^2./(pi*e*AR);

Pr = 1/2*rho*OV.^3*Cd0*S+W1^2./(0.5*rho*OV*S)*1/(pi*e*AR);

 eff_fit(c,:) = polyfit(J(:,c),eta_prop(:,c),3);
 
J_p = OV/(peaks(M)*Diam(c));

eta_curve = polyval(eff_fit(c,:),J_p);

PP(c) = max(tork(:,c).*spinny(:,c));

CP = mdot'*g/PP(c);

E = (eta_curve./CP.*sqrt(2*rho*S).*CL.^(3/2)./CD.*(W2.^(-1/2)-W1^(-1/2)));

R = eta_curve.*CL.*log(W1./W2)./(CP.*CD);

P_A = eta_curve.*PP(c);

RC = (P_A - Pr)/W1;

[beta1,RmaxTheo] = min(atan((Pr./OV)));
[PrMin,EmaxTheo] = min(Pr);
[EmaxExp,EmaxExpOP] = max(E);
[Rmax,RmaxPos] = max(R);
[RCmax,RCmaxPos] = max(RC);

store(:,g) = [(Rmax-R(RmaxTheo)) (EmaxExp-E(EmaxTheo)) RCmax];

% [bestR,RgraphNum] = min(store(1,:))
% [bestE,EgraphNum] = min(store(2,:))
% [bestRC,RCgraphNum] = max(store(3,:))

figure
plot(OV,Pr)
combo = 'Motor %d, Airframe %d, Propeller %d';
title(sprintf(combo,M,af,c))
ylabel('Power (W)')
xlabel('Airspeed (m/s)')
axis([0, max(OV), 0, max(P_A)])
hold on
plot(OV,OV*beta1)
plot(OV,P_A)
% % xline(OV(EmaxTheo),'g--');
% % xline(OV(EmaxExpOP),'g');
% % xline(OV(RmaxPos),'r');
% % xline(OV(RCmaxPos),'b');
% % xline(OV(RmaxTheo),'r--');
% 
xline((EmaxTheo),'g--');
xline((EmaxExpOP),'g');
xline((RmaxPos),'r');
xline((RCmaxPos),'b');
xline((RmaxTheo),'r--');
filename = 'm%da%dp%d.png';
saveas(gcf,sprintf(filename,M,af,c));
g = g+1;
end
end
end

[bestRAF1,RgraphNumAF1] = min(store(1,1:21));
[bestEAF1,EgraphNumAF1] = min(store(2,1:21));
[bestRCAF1,RCgraphNumAF1] = max(store(3,1:21));

[bestRAF2,RgraphNumAF2] = min(store(1,22:42));
[bestEAF2,EgraphNumAF2] = min(store(2,22:42));
[bestRCAF2,RCgraphNumAF2] = max(store(3,22:42));

