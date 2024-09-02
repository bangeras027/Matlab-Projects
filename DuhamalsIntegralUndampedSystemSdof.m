clear all
clear
clc 

% Duhamal integral for undamped system

% Input
deltaT = .1 ; % Time step for integration
tspan = 0:deltaT:15; % Time span for response 
mass = 100; % Mass of system
stiff = 100000; % Stiffness of system
omega = sqrt(stiff/mass) ; % Circular frequency

% Excitation function
for i = 1:length(tspan)
    if tspan(i) <= 10
        force(i) = - .01 * 386;
    else
        force(i) = 0;
    end
end


% Calculation of recurrent term
A(1) = 0; % Initialize
B(1) = 0;
for i = 2:length(tspan)
deltaF = force(i) - force(i-1);
deltaTime = tspan(i) - tspan(i-1);
syms tou
func1 = (force(i-1) + deltaF/deltaTime * (tou-tspan(i-1))) * cos(omega * tou) ;
func2 = (force(i-1) + deltaF/deltaTime * (tou-tspan(i-1))) * sin(omega * tou) ;
A(i) = A(i-1) + int(func1,tou,tspan(i-1),tspan(i));
B(i) = B(i-1) + int(func2,tou,tspan(i-1),tspan(i));
u(i) = (A(i)*sin(omega*tspan(i)) - B(i)*cos(omega*tspan(i)))/(omega); % Displacement 
end

figure;
plot(tspan,u);
xlabel('time');
ylabel('displacement');
axis tight
grid on





