function twoLinkRobotAdaptive
close all
global controlinput
%Set up parameters for sim
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;

% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2];
controlinput = [1;1];
% Simulation final time
tf   = 100;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e., in this sim, X0 = [e1;e2;thetahat0])
X0   = [4;10;3;2;1;1];

% Options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Integrate (you can send the paramters theta to the dynamics as seen below)
[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,opts);

% Set up desired trajectory data for plots (enter desired trajectory for your simulation)
qd = [cos(0.5*t) 2*cos(t)]';

% Parse integrated states (STATES is the same "form" as X0)
% (i.e., in this sim, STATES = [e r thetahat] over all time);
e1  = STATES(:,1:2)';
e2  = STATES(:,3:4)';
lambda1 = STATES(:,5:6)';



% Compute x from e and xd for plotting purposes
q  = e1 + qd;

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
plot(t,e2,'--','LineWidth',2)
title('Filtered Tracking error')
xlabel('time');
ylabel('Amplitude');
% 
% % Plot the adaptive estimates vs actual parameters
% figure(3)
% plot(t,repmat(theta,1,length(t)),'-','LineWidth',2)
% hold on
% ax = gca;
% ax.ColorOrderIndex = 1;
% plot(t,thetaHat,':','LineWidth',2)
% title('Adaptive Estimates & Actual Parameters')
% xlabel('time');
% ylabel('Amplitude');
% hold off

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
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Select gains for controller
K        = 5; %Enter a number
a        = 1200; %Enter a number
gamma    = 10;
alpha2   = 0.002;
beta1    = 3.6;

% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)];%Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e1          = [X(1);X(2)];
e2          = [X(3);X(4)];
lambda1     = [X(5);X(6)];

% Compute current x and xDot for convenience
q        = e1 + qd;
qDot     = e2 - a*e1 + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));


% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];

% % Compute current regression matrix
% y11      = -qdDotDot(1)+a*(qDot(1)-qdDot(1)); %Enter the expression
% y12      = -qdDotDot(2)+a*(qDot(1)-qdDot(2)); %Enter the expression
% y13      = s2*qDot(2)*qdDot(1) + s2*qDot(1)*qdDot(2) + s2*qDot(2)*qdDot(2) + a*s2*qDot(2)*e(1) + a*s2*qDot(1)*e(2) + a*s2*qDot(2)*e(2) - 2*c2*qdDotDot(1) - c2*qdDotDot(2) + 2*a*c2*(qDot(1)-qdDot(1)) + a*c2*(qDot(2)-qdDot(2)); %Enter the expression
% y14      = -qDot(1);
% y15      = 0;
% 
% y21      = 0; %Enter the expression
% y22      = -qdDotDot(1) - qdDotDot(2) + a*(qDot(1)-qdDot(1)) + a*(qDot(2)-qdDot(2)); %Enter the expression
% y23      = -s2*qDot(1)*qdDot(1) - a*s2*qDot(1)*e(1) - c2*qdDotDot(1) + a*c2*(qDot(1)-qdDot(1))   ; %Enter the expression
% y24      = 0;
% y25      = -qDot(2);
% Y        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];

% Design controller

u        = -(K+1)*e2 + (K+1)*([3;2])+lambda1; %Enter the expression
controlinput = horzcat(controlinput,u);

% Compute current closed-loop dynamics
e1Dot        = qDot - qdDot;
e2Dot        = M\(-Vm*qDot-fd*qDot+u-M*qdDotDot+M*a*e1Dot); %Enter the expression
lambda1Dot      = (K + 1)*alpha2*e2 + beta1*sign(e2);
 %Enter the expression

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [e1Dot;e2Dot;lambda1Dot];
