function [out] = Proj(theta_hat_curr, theta_hat_dot, upperlim, lowerlim,dt) %functio
[I,J] = size(theta_hat_curr); %initialy matrix definition
out = zeros(I,J);             %initial dummy matrix to be edited later
theta_hat = theta_hat_curr + dt*theta_hat_dot; %proposed update 

% check the proposed update for each index and don't allow it to change
% if it will drive the estimate out of the upper or lower limit
for i = 1:I
    for j = 1:J
        if theta_hat(i,j) >= upperlim && theta_hat_dot(i,j) > 0
            out(i,j) = 0;
        elseif theta_hat(i,j) <= lowerlim && theta_hat_dot(i,j) < 0
            out(i,j) = 0;
        else
            out(i,j) = theta_hat_dot(i,j);
        end
    end
end

end