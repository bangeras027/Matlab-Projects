clear all
clear
clc 

% Duhamal integral for damped system

data = load('el_centro.txt');

% Input
% deltaT = .1 ; % Time step for integration
tspan = data(:,1); % Time span for response 
natFreq = .1:.1:1 ; % Natural frequency cps
omega = 2*pi*natFreq ; % Circular frequency
zeta = .2 ; % Damping ratio
omegaD = omega * sqrt(1-zeta^2) ; % Damped circular frequency

% Excitation function
for k = 1:length(omega)

force = -1*data(:,2)*386; % acceleration data

% Calculation of recurrent term
Ad(1) = 0; % Initialize
Bd(1) = 0;
for i = 2:length(tspan)
deltaF = force(i) - force(i-1);
deltaTime = tspan(i) - tspan(i-1);
syms tou
func1 = (force(i-1) + deltaF/deltaTime * (tou-tspan(i-1))) * exp(zeta*omega(k)*tou) * cos(omegaD(k) * tou) ;
func2 = (force(i-1) + deltaF/deltaTime * (tou-tspan(i-1))) * exp(zeta*omega(k)*tou) *sin(omegaD(k) * tou) ;
Ad(i) = Ad(i-1) + int(func1,tou,tspan(i-1),tspan(i));
Bd(i) = Bd(i-1) + int(func2,tou,tspan(i-1),tspan(i));
u(i) = (Ad(i)*sin(omegaD(k)*tspan(i)) - Bd(i)*cos(omegaD(k)*tspan(i))) * exp(-zeta*omega(k)*tspan(i)) / (omegaD(k)); % Displacement 
end

% figure;
% plot(tspan,u);
% xlabe0l('time');
% ylabel('displacement');
% axis tight
% grid on

% Response Spectrum

specDisp(k) = max(abs(u)) ; % Spectral displacement
specVel(k) = omega(k) * specDisp(k) ; % Spectral Pseudovelocity
specAcc(k) = omega(k)^2 * specDisp(k) ; % Spectral Acceleration

clear Ad
clear Bd

end

subplot(1,3,1)
plot(natFreq,specDisp)
ylabel('SpectralDisplacement')
xlabel('Natural frequency')
axis tight
grid on

subplot(1,3,2)
plot(natFreq,specVel)
ylabel('SpectralPseudoVelocity')
xlabel('Natural frequency')
axis tight
grid on

subplot(1,3,3)
plot(natFreq,specAcc)
ylabel('SpectralAcceleration')
xlabel('Natural frequency')
axis tight
grid on
