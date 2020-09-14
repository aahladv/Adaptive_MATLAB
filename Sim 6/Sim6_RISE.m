clear
clc

%% Parameters



%% GAINS

k = 50;
kn = 10;
ks = 1;
kr = 3;

alpha = 10;
gamma1 = 5;
gamma2 = 10;
beta1  = 2;

proj1 = 4.5;   %vHat Saturation
proj2 = 6.7;    %wHat Saturation

v_initial = zeros(7,5);
w_initial = zeros(6,2);
lambda_dot_initial = (kr + 1)*alpha +[0.5;0] + beta1*[sign(0.5); sign(0)];


%% Run Simulation
[t,~,states,u,qd,f_hat,f,w_hat,v_hat] = sim('Model_Composite_2Link_NN_RISE');

%% Analysis/Plot

states = [states(:,3:4),states(:,1:2)]; % reorder states to [pos1, pos2, vel1, vel2]
error = qd-states;
rms = sqrt(mean(error.^2))*180/pi;

[dim1,dim2,dim3] = size(v_hat);
 
%Pulling Information to Plot vHat
v_hat_track = zeros(7,dim3);
for ii = 1:dim3
    index = v_hat(:,1,ii);
    v_hat_track(:,ii)=index; 
end
 
 
%Pulling Information to Plot wHat
w_hat_track = zeros(4,dim3); 
for ii = 1:dim3
    index = [w_hat(1,1,ii);
             w_hat(2,1,ii);
             w_hat(1,2,ii);
             w_hat(2,2,ii);];
    w_hat_track(:,ii) = index;
end

figure(2)
plot(t,u(:,1))
title('Link 1 Control')
xlabel('Time (s)')
ylabel('Control Torque (N-m)')
axis([0,25,-250,250])
saveas(figure(2),'Link 1 Control sm.png')

 figure(3)
 plot(t,u(:,2))
 title('Link 2 Control')
 xlabel('Time (s)')
 ylabel('Control Torque (N-m)')
 axis([0,25,-30,30])
 saveas(figure(3),'Link 2 Control sm.png')
 
figure(4)
plot(t,error(:,1)*180/pi,t,error(:,3)*180/pi)
title('Link 1 Position Errors')
legend('q1 Error','q1Dot Error')
xlabel('Time (s)')
ylabel('Link Position and Velocity Error (deg, deg/s)')
axis([0,25,-50,50])
saveas(figure(4),'Link 1 Position Errors sm.png')

 figure(5)
 plot(t,error(:,2)*180/pi,t,error(:,4)*180/pi)
 title('Link 2 Position Errors')
 legend('q2 Error','q2Dot Error')
 xlabel('Time (s)')
 ylabel('Link Position and Velocity Error (deg, deg/s)')
 axis([0,25,-50,50])
 saveas(figure(5),'Link 2 Position Errors sm.png')
 
  figure(6)
  plot(t,(f(:,1)-f_hat(:,1)))
  title('Link 1 Parameter Estimate Errors')
  xlabel('Time (s)')
  ylabel('Percent Error (%)')
  axis([0,25,-50,50])
  saveas(figure(6),'Link 1 Parameter Estimate Errors sm.png')
  
  figure(7)
  plot(t,(f(:,2)-f_hat(:,2)))
  title('Link 2 Parameter Estimate Errors')
  xlabel('Time (s)')
  ylabel('Percent Error (%)')
  axis([0,25,-50,50])
  saveas(figure(7),'Link 2 Parameter Estimate Errors sm.png')
   
figure(8)
plot(t,v_hat_track)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
title('vHat vs Time')
xlabel('Time (s)')
ylabel('vHat')
axis([0,25,-6,6])
legend('vHat1','vHat2','vHat3','vHat4','vHat5','vHat6','vHat7','Location','northeast') 
% saveas(figure(8),'wHat Graph')
 
figure(9)
plot(t,w_hat_track)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
title('wHat vs Time')
xlabel('Time (s)')
ylabel('wHat')
axis([0,25,-6,6])
legend('wHat1','wHat2','wHat3','wHat4','Location','northeast') 
% saveas(figure(8),'wHat Graph')
