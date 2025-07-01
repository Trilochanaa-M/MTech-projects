
% System parameters
k = 0.02; B = 0.5; d = 0.5; fmax = 4; w = 1; Tc = 2.4; H = 0.5915; Po = 1.26e11; La = 43655; kA = 8.7e9; Csoo = 1480; kp = 0.175; kmm = 1.478; Km = 108e-6; kc = 26e-6; kr = 0.0828; ksr = 0.0303; kb = 157.072; kt = 0.0828; X = 0.2; Fo = 2.25e-2; A = 0.203; s = 1.794; Ea = 54.63; a = 5.101e14; R = 8.314; z = 40; S = 3.9e10; c = 4.22e23; Coco = 1.4e5; fgtm = 8.3259e13; Cato = 590; ka = 1.773e20; To = 288;     
data = readtable('global_1751_2017.csv','ReadVariableNames',false);
t1 = data.Var1;
y = data.Var2;
y1 = y;

%ODE simulation
[t,pro_] = ode15s(@(t,pro)myode(t,pro),linspace(1900,2200,123),[0.01,590,1.4*10^5,540,1480,288]);

x1_pro = pro_(:,1); Cat_pro = pro_(:,2); Coc_pro = pro_(:,3); Cveg_pro = pro_(:,4); Cso_pro = pro_(:,5); T_pro = pro_(:,6);
figure
subplot(2,3,1)
plot(t,x1_pro)
title('x vs t')
xlabel('t in yrs')
ylabel('x')
subplot(2,3,2)
plot(t,Cat_pro)
title('Cat vs t')
xlabel('t in yrs')
ylabel('Cat')
subplot(2,3,3)
plot(t,Coc_pro)
title('Coc vs t')
xlabel('t in yrs')
ylabel('Coc')
subplot(2,3,4)
plot(t,Cveg_pro)
title('Cveg vs t')
xlabel('t in yrs')
ylabel('Cveg')
subplot(2,3,5)
plot(t,Cso_pro)
title('Cso vs t')
xlabel('t in yrs')
ylabel('Cso')
subplot(2,3,6)
plot(t,T_pro)
title('T vs t')
xlabel('t in yrs')
ylabel('T in K')
T_pro(end)

%ODE simulation for control
[t,pro_] = ode15s(@(t,pro)myode(t,pro),linspace(1900,2022,123),[0.01,590,1.4*10^5,540,1480,288]);
x1_pro = pro_(:,1); Cat_pro = pro_(:,2); Coc_pro = pro_(:,3); Cveg_pro = pro_(:,4); Cso_pro = pro_(:,5); T_pro = pro_(:,6);
k1=k;
T_pro(end)
%% Control 
% State  variables
syms x1 Cat Coc Cveg Cso T

% Input variables
syms et k n

% defining the equations
f = (fmax) / (1 + exp(-w*(T-Tc)));
pCO2a = (fgtm * (Cat + Cato) / ka);
tau = (1.73 * (pCO2a) ^ 0.263) + (0.0126 * (H * Po * exp(-La / (R * T))) ^ 0.503) + 0.0231;
P = n * kp * Csoo * kmm * ((pCO2a - kc) / (Km + pCO2a - kc)) * ((((15 + T) ^ 2) * (25 - T)) / 5625);
Rveg = n * kr * Cveg * kA * exp((-Ea) / (R * (T + To)));
Rso = ksr * Cso * kb * exp((-308.56) / (T + To - 227.13));
L = kt * Cveg;
Foc = Fo * X * (Cat - (z * (Cato / Coco) * Coc));
Fd = ((((1 - A) * S) / 4) * (1 + ((3 / 4) * tau)));

% defining state equations
f1 = k*x1*(1-x1)*(-B+f+(d*((2*x1)-1)));
f2 = et*(1-x1) - (P) + Rveg + (Rso)-Foc;
f3 = Foc;
f4 = P-Rveg-L;
f5 = L-Rso;
f6 = (((Fd-(s*T^4))*a)/c);

% Define the vector of variables
DX_ = [f1,f2,f3,f4,f5,f6]; % state functions
X_ = [x1,Cat,Coc,Cveg,Cso,T]; % state variables
U_ = [et,k,n]; % input variables
Y_ = [Cat T]; % measurement variables

% Compute the Jacobian matrix
A = jacobian(DX_,X_);
B = jacobian(DX_,U_);
C = jacobian(Y_,X_);
D = jacobian(Y_,U_);

% Evaluate the Jacobian at a particular point
et = e(2200); k=k1; n=1;
x1 = x1_pro(end); Cat = Cat_pro(end); Coc = Coc_pro(end); Cveg = Cveg_pro(end);Cso = Cso_pro(end);T = T_pro(end);
A = eval(subs(A)); B = eval(subs(B)); C = eval(subs(C)); D=eval(subs(D));

%% Models
%State space model
stsp1 = ss(A,B,C,D);
%Transfer function model
trfn1 = tf(stsp1);
%% function for e(t)
function y2 = e(t)
data = readtable('global_1751_2017.csv','ReadVariableNames',false);
t1 = data.Var1;
y = data.Var2;
y1 = y;
        if t <= 2014
            y2 = interp1(t1,y1,t);
        else
        %calculate emissions for t > 2014
            emax = 4.2;
            s = 30;
            y2 = y1(263)  + (((t-2014)*emax)/(t-2014+s));
        end
end

%% The function required for simulation
function dXdt = myode(t,X)
    x1 = X(1); Cat = X(2)-590; Coc = X(3)-1.4e5; Cveg = X(4)-540; Cso = X(5)-1480; T = X(6); % states % inputs
 
    % System parameters
    k = 0.02; B = 0.5; d = 0.5; fmax = 4; w = 1; Tc = 2.4; H = 0.5915; Po = 1.26e11; La = 43655; kA = 8.7e9; Csoo = 1480; kp = 0.175; kmm = 1.478; Km = 108e-6; kc = 26e-6; kr = 0.0828; ksr = 0.0303; kb = 157.072; kt = 0.0828; X = 0.2; Fo = 2.25e-2; A = 0.203; s = 1.794; Ea = 54.63; a = 5.101e14; R = 8.314; z = 40; S = 3.9e10; c = 4.22e23; Coco = 1.4e5; fgtm = 8.3259e13; Cato = 590; ka = 1.773e20; To = 288;    

    f = (fmax)/(1+exp(-w*(T-Tc)));
    pCO2a = (fgtm*(Cat+Cato)/ka);
    tau = (1.73*(pCO2a)^0.263) + (0.0126*(H*Po*exp(-La/(R*T)))^0.503) +0.0231;
    P = kp*Csoo*kmm*((pCO2a-kc)/(Km+pCO2a-kc))*((((15+T)^2)*(25-T))/5625);
    Rveg = kr*Cveg*kA*(exp((-Ea)/(R*(T+To))));
    Rso = ksr*Cso*kb*exp((-308.56)/(T+To-227.13));
    L = kt*Cveg;
    Foc = Fo*X*(Cat-(z*(Cato/Coco)*Coc));
    Fd = ((((1-A)*S)/4)*(1+((3/4)*tau)));

    dXdt(1,1) = k*x1*(1-x1)*(-B+f+(d*((2*x1)-1)));
    dXdt(2,1) = e(t)*(1-x1) - (P) +Rveg + (Rso)-Foc;
    dXdt(3,1) = Foc;
    dXdt(4,1) = P-Rveg-L;
    dXdt(5,1) = L-Rso;
    dXdt(6,1) = (((Fd-(s*T^4))*a)/c);
end