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
tf   = 10;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e1., in this sim, X0 = [e0;r0;thetahat0])
X0   = [4;10;3;2;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];

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
uf  = STATES(:,10:11)';
Ydf = STATES(:,12:13:14:15:16:17:18:19:20:21)';
lambda2 = STATES(:,22:23:24:25:26)';
mu2 = STATES(:,27:28)';
lambda1 = STATES(:,29:30)';

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
K        = 10   ; %Enter a number
a        = 20    ; %Enter a number
v = [20 5 5 40 5];
gamma    = diag(v);
beta     = 0.99;
alpha2   = 0.05;
K2       = 10;

% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression
qdDotDotDot =[0.125*sin(0.5*t); 2*sin(t)]; %Enter the expression
qdDotDotDotDot =[0.0625*cos(0.5*t); -2*cos(t)]; %Enter the expression


% Parse current states (X is the same size and "form" as X0)
% (i.e1., in this sim, X = [e1;r;thetahat])
e1        = [X(1);X(2)];
e2        = [X(3);X(4)];
thetaHat = [X(5);X(6);X(7);X(8);X(9)];
uf       = [X(10);X(11)];
Ydf      = [X(12) X(13) X(14) X(15) X(16);X(17) X(18) X(19) X(20) X(21)];
Ydf2     = [X(12);X(13);X(14);X(15);X(16);X(17);X(18);X(19);X(20);X(21)];
lambda2  = [X(22);X(23);X(24);X(25);X(26)];
mu2      = [X(27);X(28)];
lambda1  = [X(29);X(30)]; 
% Compute current x and xDot for convenience
q        = e1 + qd;
qDot     = e2 - a*e1 + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));

cd2      = cos(qd(2));
sd2      = sin(qd(2));


% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];

% Compute current regression matrix
eDot = qDot - qdDot;

b1      = qdDot(1) + a*e1(1);
b2      = qdDot(2) + a*e1(2);

b3      = qdDot(1) - qDot(1);
b4      = qdDot(2) - qDot(2);

b5      = qdDotDot(1) + a*eDot(1);
b6      = qdDotDot(2) + a*eDot(2);

% compute current regression matrix (Y2 for gradient term/ non cl/icl)
%missing%



%Y2 Matrix
y2_11   = -b5;
y2_12   = -b6;
y2_13   = -2*c2*b5 - c2*b6 + s2*qDot(2)*b1 + s2*b2*(qDot(1)+qDot(2));
y2_14   = -qDot(1);
y2_15   = 0;

y2_21   = 0;
y2_22   = -b5 - b6;
y2_23   = -c2*b5 - s2*qDot(1)*b1;
y2_24   = 0;
y2_25   = -qDot(2);

Y      = [y2_11 y2_12 y2_13 y2_14 y2_15; y2_21 y2_22 y2_23 y2_24 y2_25];



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
Yd2      = [yd11; yd12; yd13; yd14; yd15; yd21; yd22; yd23; yd24; yd25];
 %initialize the value

 %Design YdDot Matrix
 YdDot11  = qdDotDotDot(1);
 YdDot12  = qdDotDotDot(2);
 YdDot13  = 2*c2*qdDotDotDot(1) - 2*qdDotDot(1)*s2*qdDot(2) + c2*qdDotDotDot(2) - s2*qdDotDot(2)*qdDot(2) - 2*c2*qdDot(2)*qdDot(2)*qdDot(1) - 2*s2*qdDotDot(2)*qdDot(1) - 2*s2*qdDot(2)*qdDotDot(1) - c2*qdDot(2)*qdDot(2)*qdDot(2) - 2*s2*qdDotDot(2)*qdDot(2);
 YdDot14  = qdDotDotDot(1);
 YdDot15  = 0;
 
 YdDot21  = 0;
 YdDot22  = qdDotDotDot(1) + qdDotDotDot(2);
 YdDot23  = -s2*qdDot(2) + c2*qdDotDot(1) - c2*qdDot(2) + 2*s2*qdDotDot(1)*qdDot(1);
 YdDot24  = 0;
 YdDot25  = qdDotDotDot(2);

YdDot     =[YdDot11 YdDot12 YdDot13 YdDot14 YdDot15;YdDot21 YdDot22 YdDot23 YdDot24 YdDot25];


%Design YdDotDot Matrix
YdDotDot11  = qdDotDotDotDot(1);
YdDotDot12  = qdDotDotDotDot(2);
YdDotDot13  = -2*s2*qdDot(2)*qdDotDotDot(1) + 2*c2*qdDotDotDotDot(1) - 2*qdDotDotDot(1)*s2*qdDot(2) - 2*qdDotDot(1)*c2*qdDot(2)*qdDot(2) - 2*qdDotDot(1)*s2*qdDotDot(2) - s2*qdDot(2)*qdDotDotDot(2) + c2*qdDotDotDotDot(2) - c2*qdDotDot(2)*qdDot(2)*qdDot(2) - s2*qdDotDotDot(2)*qdDot(2) - s2*qdDotDot(2)*qdDotDot(2) + 2*s2*qdDot(2)*qdDot(2)*qdDot(2)*qdDot(1) - 4*c2*qdDot(2)*qdDotDot(2)*qdDot(1) - 2*c2*qdDot(2)*qdDot(1) - 2*c2*qdDot(2)*qdDotDot(2)*qdDot(1) - 2*s2*qdDotDotDot(2)*qdDot(1) - 2*s2*qdDotDot(2)*qdDotDot(1) - 2*c2*qdDot(2)*qdDot(2)*qdDotDot(1) - 2*s2*qdDotDot(2)*qdDotDot(1) - 2*s2*qdDot(2)*qdDotDotDot(1) + s2*qdDot(2)*qdDot(2)*qdDot(2)*qdDot(2) - 3*c2*qdDot(2)*qdDot(2)*qdDotDot(2) - 2*c2*qdDotDot(2)*qdDot(2)*qdDot(2) - 2*s2*qdDotDotDot(2)*qdDot(2) - 2*s2*qdDotDot(2)*qdDotDot(2); 
YdDotDot14  = qdDotDotDot(1);
YdDotDot15  = 0;

YdDotDot21  = 0;
YdDotDot22  = qdDotDotDotDot(1) + qdDotDotDotDot(2);
YdDotDot23  = -c2*qdDot(2)*qdDot(2) - s2*qdDotDot(2)*qdDotDot(1) - s2*qdDot(2)*qdDotDotDot(1) - s2*qdDot(2)*qdDotDotDot(1) + c2*qdDotDotDotDot(1) + s2*qdDot(2)*qdDot(2)*qdDot(1)*qdDot(1) - c2*qdDotDot(2)*qdDot(1)*qdDot(1) - 2*c2-qdDot(2)*qdDotDot(1)*qdDot(1) + 2*c2*qdDot(2)*qdDotDot(1)*qdDot(1) + 2*s2*qdDotDotDot(1)*qdDot(1) + 2*s2*qdDotDot(1)*qdDotDot(1);
YdDotDot24  = 0;
YdDotDot25  = qdDotDotDotDot(2);

YdDotDot    = [ YdDotDot11 YdDotDot12 YdDotDot13 YdDotDot14 YdDotDot15; YdDotDot21 YdDotDot22 YdDotDot23 YdDotDot24 YdDotDot25];
 
% Design controller
mu1 = -(K+1)*e2 + (K+1)*([3;2])+lambda1;
u        = Y*thetaHat + mu1; %Enter the expression
controlinput = horzcat(controlinput,u);

% Compute current closed-loop dynamics
e1Dot       = qDot - qdDot;
epsilon     = uf - Ydf*thetaHat - mu2;
e2Dot       = M\(-Vm*qDot-fd*qDot+u-M*qdDotDot+M*a*e1Dot); %Enter the expression


thetaHatDot = (gamma*YdDot'*e2); %Enter the expression
% lambda1Dot      = (K + 1)*alpha2*e2 + beta*sign(e2);

%Filtering
ufDot       = -beta*uf + beta*u;
YdfDot      = -beta*Ydf2 + beta*Yd2;

for i = 1:5
    YdfDot2(1,i) = YdfDot(i);
    YdfDot2(2,i) = YdfDot(5+i);
end

lambda2Dot  = gamma*(-YdDotDot'*e2 + YdDot'*alpha2*e2 + YdfDot2'*epsilon);
mu2Dot      = K2*epsilon + beta*sign(e1);

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [e1Dot;e2Dot;thetaHatDot;ufDot; YdfDot; lambda2Dot; mu2Dot; lambda1Dot];
