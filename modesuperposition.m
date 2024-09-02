function [u, udot, umax, tsol] = modesuperposition(tspan, force, modeshape, omega, zeta )
% modal superpostion method to solve for response of mdof system
% tspan = time span for required response 0:deltaT:10
% force is a matrix containing force increment acting on each dof column
% wise example see below
% force = zeros(length(tspan),2) ;  % Initializtion
% for i = 1:length(tspan) % Discretizing forcing function
%    if tspan(i) <= .25
%         force(i,1) = 1000 ;
%    elseif tspan(i) > .25 && tspan(i) <=.5
%            force(i,1) = 1000/.25 * (.5 - tspan(i));
%    else
%        force(i,1) = 0;
%    end        
% end
% force(:,2) = zeros(length(tspan),1);
% force(:,3) = zeros(length(tspan),1);
% force(:,4) = zeros(length(tspan),1);

% modeshape is a matrix containing mode shapes column wise
% omega row vector containing circular frequency
% zeta row vector containing damping in each mode

% Output discription
% u matrix containing displacement of each dof column wise
% udot matrix containing velocity of each dof column wise
% umax is row containing maximum displacement of the each dof

for k = 1:size(modeshape,2)
modalForce = 0;
for i = 1:size(force, 2)
    modalForce = modalForce + modeshape(i,k) * force(:,i);
end

forceInt = @(t) interp1(tspan, modalForce, t, 'linear') ; % Apply ODEsolver
IC = [0;0] ; % Initial condtion
odeFunc = @(t,z) [z(2); forceInt(t) - omega(k)^2 * z(1) - 2 * zeta(k) * omega(k) * z(2)] ;
[tsol, zsol] = ode45(odeFunc, tspan, IC) ;
zmax(k) = max(abs(zsol(:,1))) ; % max modal response
modeResponseU(:,k) = zsol(:,1);
modeResponseUdot(:,k) = zsol(:,2);
end

% Determine the response
% displacement
for k = 1:size(force, 2)
    for i = 1:length(modeResponseU) % Transformation from modal coordinate to real coordinate
        u(i,k) = 0;
            for j = 1:size(modeshape,2)
            u(i,k) = modeshape(k,j) * modeResponseU(i,j) + u(i,k);
            end
    end
end

% Velocity
for k = 1:size(force, 2)
    for i = 1:length(modeResponseUdot) % Transformation from modal coordinate to real coordinate
        udot(i,k) = 0;
            for j = 1:size(modeshape,2)
            udot(i,k) = modeshape(k,j) * modeResponseUdot(i,j) + udot(i,k);
            end
    end
end

% Combining modal response using SRSS
umax = sqrt(modeshape.^2 * zmax'.^2) ;
end