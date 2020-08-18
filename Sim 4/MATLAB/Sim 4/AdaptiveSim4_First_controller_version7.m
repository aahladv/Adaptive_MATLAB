function twoLinkRobotAdaptive
close all
global controlinput
%Set up parameters for sim
%just changing Regression matrix and gains

p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;


% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2];
controlinput = [1;1];

% Simulation final time
tf   = 1.5;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e1., in this sim, X0 = [e0;r0;thetahat0])
X0   = [4;10;3;2;1;1;1;1;1;1;1];

% Options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Integrate (you can send the paramters theta to the dynamics as seen below)
[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,opts);

% Set up desired trajectory data for plots (enter desired trajectory for your simulation)
qd = [cos(0.5*t) 2*cos(t)]';

% Parse integrated states (STATES is the same "form" as X0)
% (i.e1., in this sim, STATES = [e1 r thetahat] over all time);
e1  = STATES(:,1:2)';
e2  = STATES(:,3:4)';
thetaHat = STATES(:,5:9)';
mu2 = STATES(:,10:11)';

% Compute x from e1 and xd for plotting purposes
q  = e1 + qd;

% Plot the actual vs desired trajectories
figure(1)
plot(t,qd,'-','LineWidth',2)
title('Actual and Desired Trajectories');
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,q,':','LineWidth',2)
xlabel('time')
hold off

% Plot the filtered tracking error
figure(2)
plot(t,e2,'--','LineWidth',2)
title('Filtered Tracking error')
xlabel('time')

% Plot the adaptive estimates vs actual parameters
figure(3)
plot(t,repmat(theta,1,length(t)),'-','LineWidth',2)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,thetaHat,':','LineWidth',2)
xlabel('time')

title('Adaptive Estimates & Actual Parameters')
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

function [XDot] = twoLinkdynamics(t,X,theta)
global u controlinput
t
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Select gains for controller
K        = 5   ; %Enter a number
a1       = 0.9    ; %Enter a number
a2       = 0.9;
v = [100 0.5 0.5 10 0.5];
gamma    = diag(v);
beta     = 0.7;



% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression



% Parse current states (X is the same size and "form" as X0)
% (i.e1., in this sim, X = [e1;r;thetahat])
e1        = [X(1);X(2)];
e2        = [X(3);X(4)];
thetaHat = [X(5);X(6);X(7);X(8);X(9)];
mu2      = [X(10);X(11)];

% Compute current x and xDot for convenience
q        = e1 + qd;
qDot     = e2 - a1*e1 + qdDot;



% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));

cd2      = cos(qd(2));
sd2      = sin(qd(2));


% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];

% 
% eDot = qDot - qdDot;

% 
% b1      = qdDot(1) + a1*e1(1);
% b2      = qdDot(2) + a1*e1(2);
% 
% 
% b5      = qdDotDot(1) + a1*eDot(1);
% b6      = qdDotDot(2) + a1*eDot(2);

% %Y2 Matrix
% y2_11   = -b5;
% y2_12   = -b6;
% y2_13   = -2*c2*b5 - c2*b6 + s2*qDot(2)*b1 + s2*b2*(qDot(1)+qDot(2));
% y2_14   = -qDot(1);
% y2_15   = 0;
% 
% y2_21   = 0;
% y2_22   = -b5 - b6;
% y2_23   = -c2*b5 - s2*qDot(1)*b1;
% y2_24   = 0;
% y2_25   = -qDot(2);
% 
% Y      = [y2_11 y2_12 y2_13 y2_14 y2_15; y2_21 y2_22 y2_23 y2_24 y2_25];



%Design Yd Matrix
yd11     = qdDotDot(1);
yd12     = qdDotDot(2);
yd13     = 2*cd2*qdDotDot(1) + cd2*qdDotDot(2) - sd2*qdDot(2)*qdDot(1) - sd2*qdDot(1)*qdDot(2) - sd2*qdDot(2)*qdDot(2);
yd14     = qdDot(1);
yd15     = 0;

yd21     = 0;
yd22     = qdDotDot(1)+qdDotDot(2);
yd23     = cd2*qdDotDot(1) + sd2*qdDot(1)*qdDot(1);
yd24     = 0;
yd25     = qdDot(2);
Yd       = [yd11 yd12 yd13 yd14 yd15;yd21 yd22 yd23 yd24 yd25];

 

% Design controller

u        = Yd*thetaHat + mu2; %Enter the expression

controlinput = horzcat(controlinput,u);

% Compute current closed-loop dynamics
e1Dot       = qDot - qdDot;

e2Dot       = M\(-Vm*qDot-fd*qDot+u-M*qdDotDot+M*a1*e1Dot); %Enter the expression
r        = e2Dot + a2*e2;

thetaHatDot = (gamma*Yd'*e2); %Enter the expression

mu2Dot      = (K+1)*r + beta*sign(e2);

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [e1Dot;e2Dot;thetaHatDot; mu2Dot];
