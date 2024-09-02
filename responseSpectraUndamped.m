clear all
clear
clc 

% Duhamal integral for undamped system

% Input
deltaT = .1 ; % Time step for integration
tspan = 0:deltaT:15; % Time span for response 
natFreq = .1:.1:1 ; % Natural frequency cps
omega = 2*pi*natFreq ; % Circular frequency

% Excitation function

for k = 1:length(omega)

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
func1 = (force(i-1) + deltaF/deltaTime * (tou-tspan(i-1))) * cos(omega(k) * tou) ;
func2 = (force(i-1) + deltaF/deltaTime * (tou-tspan(i-1))) * sin(omega(k) * tou) ;
A(i) = A(i-1) + int(func1,tou,tspan(i-1),tspan(i));
B(i) = B(i-1) + int(func2,tou,tspan(i-1),tspan(i));
u(i) = (A(i)*sin(omega(k)*tspan(i)) - B(i)*cos(omega(k)*tspan(i)))/(omega(k)); % Displacement 
end

% figure;
% plot(tspan,u);
% xlabel('time');
% ylabel('displacement');
% axis tight
% grid on

% Response Spectrum

specDisp(k) = max(abs(u)) ; % Spectral displacement
specVel(k) = omega(k) * specDisp(k) ; % Spectral Pseudovelocity
specAcc(k) = omega(k)^2 * specDisp(k) ; % Spectral Acceleration

clear A
clear B

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





