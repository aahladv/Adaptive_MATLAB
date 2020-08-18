function twoLinkRobotAdaptive
close all
global controlinput wHat t1 t2 u 

%Set up parameters for sim
%Errors in Yd matrix are qd - q and this entire matrix has errors qd - q

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
% (i.e1., in this sim, X0 = [e0;r0;thetahat0])
X0   = [4;10;3;2;1;1];

% Options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Integrate (you can send the paramters theta to the dynamics as seen below)
[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,opts);

% Set up desired trajectory data for plots (enter desired trajectory for your simulation)
qd = [cos(0.5*t) 2*cos(t)]';

% Parse integrated states (STATES is the same "form" as X0)
% (i.e1., in this sim, STATES = [e1 r thetahat] over all time);
e  = STATES(:,1:2)';
r  = STATES(:,3:4)';
thetaHat = STATES(:,5:6)';


% Compute x from e1 and xd for plotting purposes
q  = qd - e;

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

% figure(6)
% plot(t,repmat(theta,1,length(t))-thetaHat,'-','LineWidth',2)
% title('Parameter Estimates')
% xlabel('time');
% ylabel('Amplitude');


function [XDot] = twoLinkdynamics(t,X,theta)
global u controlinput t1 t2 
t;
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
K1       = 40;
K2       = 40;

a1       = 20;    

gamma = [0.5 0; 
         0 0.5];

beta     = 0.05;

wHat1    = K1*3;
wHat2    = K2*2;
t1       = 1;
t2       = 1;



% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression



% Parse current states (X is the same size and "form" as X0)
% (i.e1., in this sim, X = [e1;r;thetahat])
e        = [X(1);X(2)];
r        = [X(3);X(4)];
thetaHat = [X(5);X(6)];


% Compute current x and xDot for convenience
q        = qd - e;
qDot     = qdDot + a1*e - r;



% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));


% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];


%Ys Matrix
ys11    =   sign(qDot(1));
ys12    =  0;

ys21    = 0;
ys22    =   sign(qDot(2));

Ys      = [ys11 ys12;ys21 ys22];

 

% Design controller




if t == t1*4*pi
   
            if abs(wHat1) < beta
                                                             %Finding the saturated values for wHat1 
            else
                wHat1 = sign(wHat1)*beta;

            end
            
            if abs(wHat2) < beta
                                                             %Finding the saturated values for  wHat2
                else 
                    wHat2 = sign(wHat2)*beta;
            end
     wHat1 = wHat1 + K1*r(1);
     wHat2 = wHat2 + K2*r(2);                                        %Updating the wHat periodically
     t1 = t1 + 1;
else
            if t == t2*2*pi
             
                   if abs(wHat2) < beta
                                                                 %Finding the saturated values for  wHat2
                    else 
                        wHat2 = sign(wHat2)*beta;
                   end
                   
                   wHat2 = wHat2 + K2*r(2);
              t2 = t2 + 1;
            else
                    wHat1 = K1*r(1);
                    wHat2 = K2*r(2);
    end
    
    
end

wHat = [wHat1; wHat2];



u        = K*r + e + wHat + Ys*thetaHat; %Enter the expression

controlinput = horzcat(controlinput,u);

% Compute current closed-loop dynamics
eDot       = r - a1*e;

rDot       = M\(M*qdDotDot + Vm*qDot + fd*qDot +[fs1 0; 0 fs2]*[sign(qDot(1)); sign(qDot(2))]+ M*a1*eDot - u); %Enter the expression


thetaHatDot = (gamma*Ys'*r); %Enter the expression

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;rDot;thetaHatDot];
