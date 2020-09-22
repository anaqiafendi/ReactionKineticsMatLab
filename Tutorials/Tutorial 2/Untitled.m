Vr = 1000; % L
V_in = 4; % L/s
Ca_in = 5; % mol/L
k = 0.01;
Ca0 = Ca_in; % mol/L
n = 1;
clf

syms Ca(t) rA t A B C D

A = -2*k;
B = -V_in/Vr;
C = V_in*Ca_in/Vr;

rA = k * Ca^2;
dCa = diff(Ca,t);
Eq1 = diff(Ca,t) == A*Ca^2 + B*Ca + C;

% Solve ODE Numerically
% Convert symbolic ODE to vector
V = odeToVectorField(Eq1);
F = matlabFunction(V,'vars',{'t','Y'});
CAsol = ode45(F,[0,1000],[Ca0]);

% Solve Symbolicly
cond = Ca(0) == Ca0;
CAsymb(t) = dsolve(Eq1,cond);

% Evaluate numerical solution on time interval
% hold on
% x = linspace(0,1000,100000);
% y = deval(CAsol,x,1);
% plot(x,y);
% hold on

% Evaluate symbolic solution on time interval
hold on
x = linspace(0,1000,100000);
FUN = matlabFunction(CAsymb);
y = feval(FUN, x);
plot(x,y);
hold on