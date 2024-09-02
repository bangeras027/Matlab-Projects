clear all
clear
clc 

% Duhamal integral for damped system

% Input
deltaT = .02 ; % Time step for integration
tspan = 0:deltaT:.1; % Time span for response 
mass = 100 ; % Mass of system
stiff = 100000 ; % Stiffness of system
omega = sqrt(stiff/mass) ; % Circular frequency
zeta = .2 ; % Damping ratio
omegaD = omega * sqrt(1-zeta^2) ; % Damped circular frequency

% Excitation function
for i = 1:length(tspan)
    if tspan(i) <= .02
        force(i) =  (120000/.02)*tspan(i);
    elseif tspan(i)>.02 && tspan(i)<=.04
        force(i) = 120000;
    elseif tspan(i)>.04 && tspan(i)<=.06
        force(i) = (120000/.02)*(.06-tspan(i));
    else
        force(i) = 0;
    end
end

% Calculation of recurrent term
Ad(1) = 0; % Initialize
Bd(1) = 0;
for i = 2:length(tspan)
deltaF = force(i) - force(i-1);
deltaTime = tspan(i) - tspan(i-1);
syms tou
func1 = (force(i-1) + deltaF/deltaTime * (tou-tspan(i-1))) * exp(zeta*omega*tou) * cos(omegaD * tou) ;
func2 = (force(i-1) + deltaF/deltaTime * (tou-tspan(i-1))) * exp(zeta*omega*tou) *sin(omegaD * tou) ;
Ad(i) = Ad(i-1) + int(func1,tou,tspan(i-1),tspan(i));
Bd(i) = Bd(i-1) + int(func2,tou,tspan(i-1),tspan(i));
u(i) = (Ad(i)*sin(omegaD*tspan(i)) - Bd(i)*cos(omegaD*tspan(i))) * exp(-zeta*omega*tspan(i)) / (mass*omegaD); % Displacement 
end

figure;
plot(tspan,u);
xlabel('time');
ylabel('displacement');
axis tight
grid on