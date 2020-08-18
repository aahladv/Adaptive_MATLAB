function twoLinkRobotCLandICLstructure
close all
global ui u controlinput pos timeprev qDotprev Y1hist Y3prev sumY_ICL script_ui script_yi Y4 Y3 uprev

% setup and initialization
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;

pos         = 0;
timeprev    = 0;
qDotprev    = 0;
Y1hist      = [0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
ui          = 0;
       
%Stacked Paraeter vector
theta   =[p1;p2;p3;f1;f2];
controlinput = [1;1];

%simulation final time
tf  = 0.2;

%Inital condition vector
X0   = [4;10;3;2;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];

% options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% integrate (you can send the paramters theta to the dynamics as seen below)
[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,opts);

% Parse integrated states (STATES is the same "form" as X0)
% (i.e., in this sim, STATES = [e r thetahat] over all time);
e  = STATES(:,1:2)';
r  = STATES(:,3:4)';
thetaHat = STATES(:,5:9)';
Y4 = STATES(:,10:19)';
Tau = STATES(:,20:21)';
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
plot(t,repmat(theta,1,length(t)),'-r','LineWidth',2)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,thetaHat,':k','LineWidth',2)
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
global ui u controlinput pos timeprev qDotprev Y1hist Y3prev sumY_ICL script_ui script_yi Y4 Y3 uprev
t
% compute desired trajectory and needed derivatives (i.e., qd, qdDot, qdDotDot)
qd          = [cos(0.5*t) 2*cos(t)]';    
qdDot       = [-0.5*sin(0.5*t) -2*sin(t)]'; 
qdDotDot    = [-0.25*cos(0.5*t) -2*cos(t)]';

% select gains (i.e., K, Kcl/Kicl, alpha, Gamma)
K       = 10;
Kicl     = 0.0000001;
gamma   = [20 0 0 0 0;0 2 0 0 0;0 0 3 0 0;0 0 0 50 0;0 0 0 0 3]; 
a       = .3;

% parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])

e               = [X(1);X(2)];
r               = [X(3);X(4)];
thetaHat        = [X(5);X(6);X(7);X(8);X(9)];
Y4              = [X(10) X(11) X(12) X(13) X(14);X(15) X(16) X(17) X(18) X(19)];
Tau             = [X(20);X(21)];
Y42             = [X(10); X(11); X(12); X(13); X(14);X(15); X(16); X(17); X(18); X(19)];

% compute current q and qDot for convenience
q       = e + qd; 
qDot     = r - a*e + qdDot;

% compute cos(x2) and sin(x2) for convenience

c2       = cos(q(2));
s2       = sin(q(2));

cdev2   = -sin(q(2))*qDot(2); %deriving c2 and s2 for mDot term
sdev2   =  cos(q(2))*qDot(2);

cd2     = cos(qd(2));
sd2     = sin(qd(1));


%p and f values
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;

% compute current matrices for the dynamics (i.e., M, Vm, ff)


M       = [ p1+2*p3*c2 p2+p3*c2; p2+p3*c2 p2];
MDot    = [ p1+2*p3*cdev2 p2+p3*cdev2; p2+p3*cdev2 p2];
Vm      = [-p3*s2*qDot(2) -p3*s2*(qDot(1)+qDot(2)); p3*s2*qDot(1) 0];
fd      = [ f1 0; 0 f2];

eDot = qdDot - qDot;

b1      = qdDot(1) + a*e(1);
b2      = qdDot(2) + a*e(2);

b3      = qdDot(1) - qDot(1);
b4      = qdDot(2) - qDot(2);

b5      = qdDotDot(1) + a*eDot(1);
b6      = qdDotDot(2) + a*eDot(2);

% compute current regression matrix (Y2 for gradient term/ non cl/icl)
%missing%
%Y1 matrix
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
Y1        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];

%Y2 matrix

y2_11   = b5;
y2_12   = b6;
y2_13   = 2*c2*b5+c2*b6-s2*qDot(2)*b1 -s2*b2*(qDot(1)+qDot(2));
y2_14   = qDot(1);
y2_15   = 0;

y2_21   = 0;
y2_22   = b5 + b6;
y2_23   = c2*b5 + s2*qDot(1)*b1;
y2_24   = 0;
y2_25   = qDot(2);

Y2      = [y2_11 y2_12 y2_13 y2_14 y2_15; y2_21 y2_22 y2_23 y2_24 y2_25]



%Y4matrix

y4_11 = -qDot(1);
y4_12 = -qDot(2);
y4_13 = -2*cdev2*qDot(1) - cdev2*qDot(2) -s2*qDot(2)*qDot(1)- s2*qDot(2)*(qDot(1)+qDot(2));
y4_14 = qDot(1);
y4_15 = 0;

y4_21 = 0;
y4_22 = -qDot(1) - qDot(2);
y4_23 = -cdev2*qDot(1) + s2*qDot(1)*qDot(1);
y4_24 = 0;
y4_25 = qDot(2);

Y4      = [y4_11 y4_12 y4_13 y4_14 y4_15 ;y4_21 y4_22 y4_23 y4_24 y4_25]; 

% set lambda for history stack min eigenvalue condition
lambda      = 0;


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

u = Y2 * theta + e + K*r;
controlinput = horzcat(controlinput,u);

% determine which update law to use (i.e., thetaHatDot)
if pos == 1
    %(after eigenvalues condition met for Y1hist)
    
    % for CL:
    % compute new cl_sumY (function of: Yi,ui,thetaHati), 
    % assuming new data in history stack
    sumY_ICL = script_yi' * (script_ui - script_ui*thetaHat);
    
    % for ICL:
    % compute new icl_sumY(script) (function of: Yi(script),ui,x(ti),x(ti-dt),thetaHati),
    % assuming new data in history stack

    % compute thetaHatDot
    
    thetaHatDot = gamma*Y2'*r + Kicl*gamma*sumY_ICL;
else
    % (before eigenvalues condition met for Y1hist)

    % compute thetaHatDot (no cl/icl term)
    
    thetaHatDot = gamma*Y2'*r;
    
    % for CL:
    % compute cl_sumY (function of: Yi,ui, thetaHat) for use later
    %missing%
    size(sumY_ICL)
    size(gamma)
    size(Y2)
   thetaHatDot = gamma*Y2'*r + Kicl*gamma*sumY_ICL;
    % for ICL:
    % compute icl_sumY(script) (function of: Yi(script),ui,x(ti),x(ti-deltat),thetaHati) for use later
    %missing%
    
    Y3 = Y3 - Y3prev;
    
end 

% compute current closed-loop errors for integration(i.e., eDot, rDot)

eDot = qDot - qdDot;
rDot = M\(-Vm*qDot-fd*qDot+u-M*qdDotDot+M*a*eDot);

Y3          = Y3 - Y3prev;
script_yi   = Y3 + Y42;
script_ui   = u - uprev;

% Y4          = -MDot*qDot + Vm*qDot + fd*qDot;


% update "previous" variables for next loop
%(i.e., qDotprev, timeprev, Y1iprev, Uiprev)

qDotprev    = qDot;
timeprev    = t;
ui          = u;


% Stacked dynamics vector (XDot is the same size and "form" as X)

XDot        = [eDot;rDot;thetaHatDot;Y42;ui];
end

