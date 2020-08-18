function twoLinkRobotCLandICLstructure
close all
global controlinput
% setup and initialization
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;
       
%Stacked Paraeter vector
theta   =[p1;p2;p3;f1;f2];
controlinput = [1;1];

%simulation final time
tf  = 0.15;

%Inital condition vector
X0   = [4;10;3;2;1;1;1;1;1];

% options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% integrate (you can send the paramters theta to the dynamics as seen below)
[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,opts);

% Parse integrated states (STATES is the same "form" as X0)
% (i.e., in this sim, STATES = [e r thetahat] over all time);
e  = STATES(:,1:2)';
r  = STATES(:,3:4)';
thetaHat = STATES(:,5:9)';
% post computation for figures
qd          = [cos(0.5*t) 2*cos(t)]';
q = e + qd;

% plot figures
% Plot the actual vs desired trajectories
figure(1)
plot(t,qd,'-','LineWidth',2)
title('Actual and Desired Trajectories');
xlabel('time');
ylabel('Amplitude');
legend();
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,q,':','LineWidth',2)
hold off

% Plot the filtered tracking error
figure(2)
plot(t,r,'--','LineWidth',2)
title('Filtered Tracking error')
xlabel('time');
ylabel('Amplitude');

% Plot the adaptive estimates vs actual parameters
figure(3)
plot(t,repmat(theta,1,length(t)),'-','LineWidth',2)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,thetaHat,':','LineWidth',2)
title('Adaptive Estimates & Actual Parameters')
xlabel('time');
ylabel('Amplitude');
hold off

figure(4)
time = 0:tf/length(controlinput):tf;
time(end) = [];
plot(time,controlinput(1,:))
title('Control input for Link 1');
xlabel('Time(s)')
ylabel('Torque(Nm)')

figure(5)
time = 0:tf/length(controlinput):tf;
time(end) = [];
plot(time,controlinput(2,:))
title('Control input for Link 2');
xlabel('Time(s)')
ylabel('Torque(Nm)')


end

function [XDot] = twoLinkdynamics(t,X,theta)
% initialize global variables for use in plot later.
global u controlinput
t
% compute desired trajectory and needed derivatives (i.e., qd, qdDot, qdDotDot)
qd          = [cos(0.5*t) 2*cos(t)]';    
qdDot       = [-0.5*sin(0.5*t) -2*sin(t)]'; 
qdDotDot    = [-0.25*cos(0.5*t) -2*cos(t)]';

% select gains (i.e., K, Kcl/Kicl, alpha, Gamma)
K       = 10;
Kcl     = 12;
gamma   = 12;
a       = 13;

% parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])

e               = [X(1);X(2)];
r               = [X(3);X(4)];
thetaHat        = [X(5);X(6);X(7);X(8);X(9)];

% compute current q and qDot for convenience
q       = e + qd; 
qDot     = r - a*e + qdDot;
% compute cos(x2) and sin(x2) for convenience

c2       = cos(q(2));
s2       = sin(q(2));

cd2     = cos(qd(2));
sd2     = sin(qd(1));


%p and f values
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;

% compute current matrices for the dynamics (i.e., M, Vm, ff)


M   = [ p1+2*p3*c2 p2+p3*c2; p2+p3*c2 p2];
Vm  = [-p3*s2*qDot(2) -p3*s2*(qDot(1)+qDot(2)); p3*s2*qDot(1) 0];
fd  = [ f1 0; 0 f2];


% compute current regression matrix (Y2 for gradient term/ non cl/icl)
%missing%
%Y2 matrix
y11      = -qdDotDot(1)+a*(qDot(1)-qdDot(1)); %Enter the expression
y12      = -qdDotDot(2)+a*(qDot(2)-qdDot(2)); %Enter the expression
y13      = s2*qDot(2)*qdDot(1) + s2*qDot(1)*qdDot(2) + s2*qDot(2)*qdDot(2) + a*s2*qDot(2)*e(1) + a*s2*qDot(1)*e(2) + a*s2*qDot(2)*e(2) - 2*c2*qdDotDot(1) - c2*qdDotDot(2) + 2*a*c2*(qDot(1)-qdDot(1)) + a*c2*(qDot(2)-qdDot(2)); %Enter the expression
y14      = -qDot(1);
y15      = 0;

y21      = 0; %Enter the expression
y22      = -qdDotDot(1) - qdDotDot(2) + a*(qDot(1)-qdDot(1)) + a*(qDot(2)-qdDot(2)); %Enter the expression
y23      = -s2*qDot(1)*qdDot(1) - a*s2*qDot(1)*e(1) - c2*qdDotDot(1) + a*c2*(qDot(1)-qdDot(1))   ; %Enter the expression
y24      = 0;
y25      = -qDot(2);
Y2        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];

%Y2d matrix 
yd11    = qdDotDot(1);
yd12    = qdDotDot(2);
yd13    = 2*cd2*qdDotDot(1) + cd2*qdDotDot(2) - sd2*qdDot(2)*qdDot(1) - s2*qdDot(1)*qdDot(2) - s2*qdDot(2)*qdDot(2);
yd14    = qdDot(1);
yd15    = 0;

yd21    = 0;
yd22    = qdDotDot(1) + qdDotDot(2);
yd23    = cd2*qdDotDot(1) + s2*qdDot(2)*qd(1);
yd24    = 0;
yd25    = qdDot(2);

Y2d = [yd11 yd12 yd13 yd14 yd15;yd21 yd22 yd23 yd24 yd25];



% set lambda for history stack min eigenvalue condition
lambda      = 6;
pos         = 0;
timeprev    = 0;
qDotprev    = 0;
Y1hist      = [0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
ui          = 0;

% determine if the history stack meets the min eigenvalue condition
if pos == 0 %if min eigenvalue of history stack less than lambda

    % compute qDotDot for Y1i
    dt        = t - timeprev; % determine delta time.
    dq        = qDot - qDotprev; % determine delta qDot.
    qDotDot   = dq/dt; % determine qDotDot.

    %if first loop (cannot determine dq or dt), initiallize qDotDot.
    if dt == 0
        qDotDot   = [0;0];
    end

    %determine Y1i (Y1 for this loop)
    %missing%
    Y1i     = [0 0 0 0 0;0 0 0 0 0];
    
    % determine new history stack (summation)
    Y1hist    = Y1hist + Y1i'*Y1i;
    
    % optional: save the Y1i data in an array instead, then add/remove data
    % from history stack to maximize min eigenvalue of Y1hist.

    % determine min eigenvalue of history stack
    Y1eigenVal = eig(Y1hist);

    % check min eigenvalue condition
    Y1eigenValMin = min(Y1eigenVal);
    if (Y1eigenValMin > lambda)
       pos = 1;
    else
       pos = 0;
    end
end

% design controller (i.e., u)

u = Y2 * theta;
controlinput = horzcat(controlinput,u);

% determine which update law to use (i.e., thetaHatDot)
if pos == 1
    %(after eigenvalues condition met for Y1hist)
    
    % for CL:
    % compute new cl_sumY (function of: Yi,ui,thetaHati), 
    % assuming new data in history stack
    c1_sumY = Y1i' * (ui - Y1i*thetaHat);
    
    % for ICL:
    % compute new icl_sumY(script) (function of: Yi(script),ui,x(ti),x(ti-dt),thetaHati),
    % assuming new data in history stack

    % compute thetaHatDot
    
    thetaHatDot = gamma*Y2'*r + Kcl*gamma*c1_sumY;
else
    % (before eigenvalues condition met for Y1hist)

    % compute thetaHatDot (no cl/icl term)
    
    thetaHatDot = gamma*Y2'*r;
    
    % for CL:
    % compute cl_sumY (function of: Yi,ui, thetaHat) for use later
    %missing%
    
    c1_sumY = Y1i' * (ui - Y1i*thetaHat);
    % for ICL:
    % compute icl_sumY(script) (function of: Yi(script),ui,x(ti),x(ti-deltat),thetaHati) for use later
    %missing%
end 

% compute current closed-loop errors for integration(i.e., eDot, rDot)

eDot = qDot - qdDot;
rDot = M\(Vm*(qdDot + a*e) +fd*qDot +a*eDot - u - Vm*r + M*qdDotDot);
% update "previous" variables for next loop
%(i.e., qDotprev, timeprev, Y1iprev, Uiprev)

qDotprev    = qDotprev  +   1;
timeprev    = timeprev  +   1;
Y1i         = Y1i       +   1;
ui          = ui        +   1;

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;rDot;thetaHatDot];
end

