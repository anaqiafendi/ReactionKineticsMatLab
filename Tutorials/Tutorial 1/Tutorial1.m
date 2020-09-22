% Solve the ODE using matlab

% https://blogs.mathworks.com/loren/2013/06/10/from-symbolic-differential-equations-to-their-numeric-solution/

Vr = 1000; % L
V0 = 5; % L/s
k = 0.01;
CA0 = 6; % mol/L
n = 1;
clf
% CA0 = 8;

for n = [-1,0,0.5,1,1.5,2]
    % Establish ODE Symbolically
    syms CA(t) rA t
    rA = -k * CA^n;
    dCA = diff(CA,t);
    Eq1 = Vr*diff(CA,t) == CA0*V0 - CA*V0 + Vr*rA;
    
    % Convert symbolic ODE to vector
    V = odeToVectorField(Eq1);
    F = matlabFunction(V,'vars',{'t','Y'});
    CAsol = ode45(F,[0,1000],[CA0]);
    
    % Evaluate solution on time interval
    hold on
    x = linspace(0,1000,10000);
    y = deval(CAsol,x,1);
    plot(x,y);
    hold on
end

legend('n = -1','n = 0','n = 0.5','n = 1','n = 1.5','n = 2');
ylabel('Concentration (mol/L)')
xlabel('Time (s)')
title('Concentration over time of A across different n values')