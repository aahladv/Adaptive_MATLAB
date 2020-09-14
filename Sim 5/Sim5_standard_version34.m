function twoLinkRobotAdaptive
close all
global controlinput u

%Set up parameters for sim
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;
fs1      = 1.2;
fs2      = 0.4;


% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2;fs1;fs2];
controlinput = [1;1]; 

% Simulation final time
tf   = 60;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e., in this sim, X0 = [e0;r0;thetahat0])
X0   = [4;10;3;2;1;1;1;1;1;1;1];


% Options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Integrate (you can send the paramters theta to the dynamics as seen below)
[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,opts);

% Set up desired trajectory data for plots (enter desired trajectory for your simulation)
qd = [cos(0.5*t) 2*cos(t)]';
 
% Parse integrated states (STATES is the same "form" as X0)
% (i.e., in this sim, STATES = [e r thetahat] over all time);
e  = STATES(:,1:2)';
r  = STATES(:,3:4)';
thetaHat = STATES(:,5:6:7:8:9:10:11)';

% Compute x from e and xd for plotting purposes
q  = e + qd;



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
plot(t,r,'--','LineWidth',2)
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

figure(6)
plot(t,repmat(theta,1,length(t))-thetaHat,'-','LineWidth',2)
title('Parameter Estimates')
xlabel('time');
ylabel('Amplitude');



function [XDot] = twoLinkdynamics(t,X,theta)
global u controlinput
% Parse parameter vector
p1  = theta(1);
p2  = theta(2);
p3  = theta(3);
f1  = theta(4);
f2  = theta(5);
fs1 = theta(6);
fs2 = theta(7);

% Select gains for controller
K        = 40;  
a        = 0.5; 
ddd      = [ 40 10.50 0 40 5 0.001 0.0001];

gamma    = diag(ddd);


%Alternative gains % ddd      = [ 39 5 15 50 5 0.001 0.0001];


% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)];  
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)];  

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e        = [X(1);X(2)];
r        = [X(3);X(4)];
thetaHat = [X(5);X(6);X(7);X(8);X(9);X(10);X(11)];


% Compute current x and xDot for convenience
q        = e + qd;
qDot     = r - a*e + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));




% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];

% Compute current regression matrix
y11      = -qdDotDot(1)+a*(qDot(1)-qdDot(1)); 
y12      = -qdDotDot(2)+a*(qDot(1)-qdDot(2)); 
y13      = s2*qDot(2)*qdDot(1) + s2*qDot(1)*qdDot(2) + s2*qDot(2)*qdDot(2) + a*s2*qDot(2)*e(1) + a*s2*qDot(1)*e(2) + a*s2*qDot(2)*e(2) - 2*c2*qdDotDot(1) - c2*qdDotDot(2) + 2*a*c2*(qDot(1)-qdDot(1)) + a*c2*(qDot(2)-qdDot(2));  
y14      = -qDot(1);
y15      = 0;
y16      = sign(qDot(1));
y17      = 0;   

y21      = 0; 
y22      = -qdDotDot(1) - qdDotDot(2) + a*(qDot(1)-qdDot(1)) + a*(qDot(2)-qdDot(2));  
y23      = -s2*qDot(1)*qdDot(1) - a*s2*qDot(1)*e(1) - c2*qdDotDot(1) + a*c2*(qDot(1)-qdDot(1))   ;  
y24      = 0;
y25      = -qDot(2);
y26      = 0;
y27      = sign(qDot(2));
Y        = [y11 y12 y13 y14 y15 y16 y17;y21 y22 y23 y24 y25 y26 y27];





% Design controller

u        = -K*r-Y*thetaHat-e; 
controlinput = horzcat(controlinput,u);

% Compute current closed-loop dynamics
eDot        = qDot - qdDot;

% size(thetaHat)
rDot        = M\(-Vm*qDot-fd*qDot+u-M*qdDotDot - [fs1 0; 0 fs2]*[sign(qDot(1)); sign(qDot(2))]+M*a*eDot);  

thetaHatDot = (gamma*Y'*r);  


% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;rDot;thetaHatDot];
