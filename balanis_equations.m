%% Balanis Equations for a Pyramidal Horn Antenna.
% Tomás Ruiz Aybar, AATM.

%% Ka-band Pyramidal Horn:

f = 32.5e9;
lambda = 3e8 / f;
a1 = 68.5e-3;
b1 = 56.5e-3;
D = sqrt(a1^2+b1^2);
p_e = 150e-3;
p_h = p_e;
b = 3.556e-3;
a = 7.112e-3;

% Angular variables
dthe = 1;    dphi = 2; % angular increments 
phi = (0:dphi:360)*(pi/180); % from 0° to 360°
theta = (0:dthe:90)*(pi/180); % from 0° to 90°
[THE, PHI] = meshgrid(theta, phi);

% Axial lenghts, use 13-47a/b
rho_e = b1*sqrt((p_e/(b1-b))^2 + 1/4);
rho_h = a1*sqrt((p_h/(a1-a))^2 + 1/4);

%phi_e = rad2deg(asin((b1/2)/rho_e));
%phi_h = rad2deg(asin((a1/2)/rho_h));
%rho_1 = rho_e*cosd(phi_e);
rho_1 = sqrt(rho_e^2 - (b1^2/4));
%rho_2 = rho_h*cosd(phi_h);
rho_2 = sqrt(rho_h^2 - (a1^2/4));

E0 = 1;
r = 2*(2*D^2)/lambda; % Verify FF conditions
k = (2*pi)/lambda;

% This section has a high cost-computational cost

% eq.[13-25, 13-26]
k_primeX = k*sin(THE).*cos(PHI) + (pi/a1);
t2_prime = sqrt(1/(pi*k*rho_2))*((k*a1)/2 - k_primeX*rho_2);
t1_prime = sqrt(1/(pi*k*rho_2))*(-(k*a1)/2 - k_primeX*rho_2);
Ct2_prime = fresnelc(t2_prime);
Ct1_prime = fresnelc(t1_prime);
St2_prime = fresnels(t2_prime);
St1_prime = fresnels(t1_prime);

k_2primeX = k*sin(THE).*cos(PHI) - (pi/a1);
t2_2prime = sqrt(1/(pi*k*rho_2))*((k*a1)/2 - k_2primeX*rho_2);
t1_2prime = sqrt(1/(pi*k*rho_2))*(-(k*a1)/2 - k_2primeX*rho_2);
Ct2_2prime = fresnelc(t2_2prime);
Ct1_2prime = fresnelc(t1_2prime);
St2_2prime = fresnels(t2_2prime);
St1_2prime = fresnels(t1_2prime);

% eq.[13-5, 13-8]
k_Y = k*sin(THE).*sin(PHI);
t1 = sqrt(1/(pi*k*rho_1))*(-(k*b1)/2 - k_Y*rho_1);
t2 = sqrt(1/(pi*k*rho_1))*((k*b1)/2 - k_Y*rho_1);
Ct2 = fresnelc(t2);
Ct1 = fresnelc(t1);
St2 = fresnels(t2);
St1 = fresnels(t1);

% eq.[13-44]
I1 = 0.5*sqrt((pi*rho_2)/k).*(exp(j.*(k_primeX.^2*rho_2)/(2*k)).*((Ct2_prime-Ct1_prime)-j.*(St2_prime-St1_prime))...
    +exp(j.*(k_2primeX.^2*rho_2)/(2*k)).*((Ct2_2prime-Ct1_2prime)-j.*(St2_2prime-St1_2prime)));

% eq.[13-45]
I2 = sqrt((pi*rho_1)/k).*exp(j.*(k_Y.^2*rho_1)/(2*k)).*((Ct2-Ct1)-j*(St2-St1));

% E field components 
Ethe = j*(k*E0*exp(-j*k*r))/(4*pi*r).*(sin(PHI).*(1+cos(THE)).*I1.*I2);
Ethe = Ethe./max(max(abs(Ethe)));

Ephi = j*(k*E0*exp(-j*k*r))/(4*pi*r).*(cos(PHI).*(cos(THE)+1).*I1.*I2);
Ephi = Ephi./max(max(abs(Ephi)));

% Plotting of far-field patterns

% Ludwig Co-Cx components
phi0 = 90*(pi/180);
Eco = Ethe.*cos(PHI-phi0) - Ephi.*sin(PHI-phi0);
Ecx = Ethe.*sin(PHI-phi0) + Ephi.*cos(PHI-phi0);

% Angular plotting variables
vthe = -90:dthe:90;
%vthe = [-fliplr(THE(1,2:end)), THE(1,:)];
vphi = 0:dphi:360-dphi;

% Find phi-indeces of E/H main planes
phi000 = 1;
phi090 = find(abs(vphi-90)<1e-3);
phi180 = find(abs(vphi-180)<1e-3);
phi270 = find(abs(vphi-270)<1e-3);
%phi090 = find(abs(PHI(:,1)-90)<1e-3);
%phi180 = find(abs(PHI(:,1)-180)<1e-3);
%phi270 = find(abs(PHI(:,1)-270)<1e-3);

% Co-polar and cx-polar phi=0 (H-plane)
Eco00 = [fliplr(Eco(phi180, 2:end)) Eco(phi000, :)];
Ecx00_t = [fliplr(Ecx(phi180, 2:end)) Ecx(phi000, :)];

% Co-polar and cx-polar phi=90 (E-plane)
Eco90 = [fliplr(Eco(phi270, 2:end)) Eco(phi090, :)];
Ecx90_t = [fliplr(Ecx(phi270, 2:end)) Ecx(phi090, :)];

% Directivity
%D = Directivity(Ethe, Ephi, vthe, vphi);
%D0 = 10*log10(max(max(D)));
u = (1/sqrt(2))*((sqrt(lambda*rho_2))/a1 + (a1/sqrt(lambda*rho_2)));
v = (1/sqrt(2))*((sqrt(lambda*rho_2))/a1 - (a1/sqrt(lambda*rho_2)));
Cu = fresnelc(u);
Cv = fresnelc(v);
Su = fresnels(u);
Sv = fresnels(v);
aux = b1/(sqrt(2*lambda*rho_1));
C2 = fresnelc(aux)^2;
S2 = fresnels(aux)^2;
Dp = (8*pi*rho_1*rho_2)/(a1*b1)*((Cu-Cv)^2 + (Su - Sv)^2)*(C2 + S2);
D0 = 10*log10(max(max(Dp)));

% COMPARISON WITH JAVIER THEORICAL RESULTS:

Balanis = 'Ka_Band_SGH_Balanis_40_00.mat';

load (Balanis)

% Ludwig Co-Cx components
phi0 = 90;
Eco = Ethe.*cosd(PHI-phi0) - Ephi.*sind(PHI-phi0);
Ecx = Ethe.*sind(PHI-phi0) + Ephi.*cosd(PHI-phi0);

% Angular plotting variables
vthe_javier = [-fliplr(THE(1,2:end)), THE(1,:)];

% Find phi-indeces of E/H main planes
phi000 = 1;
phi090 = find(abs(PHI(:,1)-90)<1e-3);
phi180 = find(abs(PHI(:,1)-180)<1e-3);
phi270 = find(abs(PHI(:,1)-270)<1e-3);

% Co-polar main planes
Eco00_javier = [fliplr(Eco(phi180, 2:end)) Eco(phi000, :)];
Eco90_javier = [fliplr(Eco(phi270, 2:end)) Eco(phi090, :)];


figure(); 
hold on; 
plot(vthe, D0+20*log10(abs(Eco00)), LineWidth=2)
plot(vthe, D0+20*log10(abs(Eco90)), LineWidth=2)
%plot(vthe, D0+20*log10(abs(Ecx00_t)), LineWidth=2)
%plot(vthe, D0+20*log10(abs(Ecx90_t)), LineWidth=2)
%plot(vthe_javier, 10*log10(abs(Eco00_javier)), 'b:' ,LineWidth=2)
%plot(vthe_javier, 10*log10(abs(Eco90_javier)), 'r:',LineWidth=2)
grid on; 
xlim([-90 90])
ylabel('Directivity [dB]', 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('\theta [º]', 'FontName', 'Times New Roman', 'FontSize', 14)
%legend('Copolar H-Plane, \phi = 0º', 'Copolar E-Plane, \phi = 90º',...
    %'Copolar H-Plane, \phi = 0º', 'Copolar E-Plane, \phi = 90º','Javier H-Plane', 'Javier E-Plane',...
    %'FontName', 'Times New Roman', 'FontSize', 14)
legend('Copolar H-Plane, \phi = 0º', 'Copolar E-Plane, \phi = 90º',...
    'FontName', 'Times New Roman', 'FontSize', 14)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

GHz32TheoretycalCoHPlane= D0+20*log10(abs(Eco00));
save('32_5GHz_theoretycal_Co_HPlane.mat', 'GHz32TheoretycalCoHPlane');

%% C-band Pyramidal Horn:

f = 4.9e9;
lambda = 3e8 / f;
a1 = 216e-3;
b1 = 160e-3;
D = sqrt(a1^2+b1^2);
p_e = 240e-3;
p_h = p_e;
b = 22.2e-3;
a = 47.55e-3;

% Angular variables
dthe = 1;    dphi = 2; % angular increments 
phi = (0:dphi:360)*(pi/180); % from 0° to 360°
theta = (0:dthe:90)*(pi/180); % from 0° to 90°
[THE, PHI] = meshgrid(theta, phi);

% Axial lenghts, use 13-47a/b
rho_e = b1*sqrt((p_e/(b1-b))^2 + 1/4);
rho_h = a1*sqrt((p_h/(a1-a))^2 + 1/4);

%phi_e = rad2deg(asin((b1/2)/rho_e));
%phi_h = rad2deg(asin((a1/2)/rho_h));
%rho_1 = rho_e*cosd(phi_e);
rho_1 = sqrt(rho_e^2 - (b1^2/4));
%rho_2 = rho_h*cosd(phi_h);
rho_2 = sqrt(rho_h^2 - (a1^2/4));

E0 = 1;
r = 2*(2*D^2)/lambda; % Verify FF conditions
k = (2*pi)/lambda;

% This section has a high cost-computational cost

% eq.[13-25, 13-26]
k_primeX = k*sin(THE).*cos(PHI) + (pi/a1);
t2_prime = sqrt(1/(pi*k*rho_2))*((k*a1)/2 - k_primeX*rho_2);
t1_prime = sqrt(1/(pi*k*rho_2))*(-(k*a1)/2 - k_primeX*rho_2);
Ct2_prime = fresnelc(t2_prime);
Ct1_prime = fresnelc(t1_prime);
St2_prime = fresnels(t2_prime);
St1_prime = fresnels(t1_prime);

k_2primeX = k*sin(THE).*cos(PHI) - (pi/a1);
t2_2prime = sqrt(1/(pi*k*rho_2))*((k*a1)/2 - k_2primeX*rho_2);
t1_2prime = sqrt(1/(pi*k*rho_2))*(-(k*a1)/2 - k_2primeX*rho_2);
Ct2_2prime = fresnelc(t2_2prime);
Ct1_2prime = fresnelc(t1_2prime);
St2_2prime = fresnels(t2_2prime);
St1_2prime = fresnels(t1_2prime);

% eq.[13-5, 13-8]
k_Y = k*sin(THE).*sin(PHI);
t1 = sqrt(1/(pi*k*rho_1))*(-(k*b1)/2 - k_Y*rho_1);
t2 = sqrt(1/(pi*k*rho_1))*((k*b1)/2 - k_Y*rho_1);
Ct2 = fresnelc(t2);
Ct1 = fresnelc(t1);
St2 = fresnels(t2);
St1 = fresnels(t1);

% eq.[13-44]
I1 = 0.5*sqrt((pi*rho_2)/k).*(exp(j.*(k_primeX.^2*rho_2)/(2*k)).*((Ct2_prime-Ct1_prime)-j.*(St2_prime-St1_prime))...
    +exp(j.*(k_2primeX.^2*rho_2)/(2*k)).*((Ct2_2prime-Ct1_2prime)-j.*(St2_2prime-St1_2prime)));

% eq.[13-45]
I2 = sqrt((pi*rho_1)/k).*exp(j.*(k_Y.^2*rho_1)/(2*k)).*((Ct2-Ct1)-j*(St2-St1));

% E field components 
Ethe = j*(k*E0*exp(-j*k*r))/(4*pi*r).*(sin(PHI).*(1+cos(THE)).*I1.*I2);
Ethe = Ethe./max(max(abs(Ethe)));

Ephi = j*(k*E0*exp(-j*k*r))/(4*pi*r).*(cos(PHI).*(cos(THE)+1).*I1.*I2);
Ephi = Ephi./max(max(abs(Ephi)));

% Plotting of far-field patterns

% Ludwig Co-Cx components
phi0 = 90*(pi/180);
Eco = Ethe.*cos(PHI-phi0) - Ephi.*sin(PHI-phi0);
Ecx = Ethe.*sin(PHI-phi0) + Ephi.*cos(PHI-phi0);

% Angular plotting variables
vthe = -90:dthe:90;
%vthe = [-fliplr(THE(1,2:end)), THE(1,:)];
vphi = 0:dphi:360-dphi;

% Find phi-indeces of E/H main planes
phi000 = 1;
phi090 = find(abs(vphi-90)<1e-3);
phi180 = find(abs(vphi-180)<1e-3);
phi270 = find(abs(vphi-270)<1e-3);
%phi090 = find(abs(PHI(:,1)-90)<1e-3);
%phi180 = find(abs(PHI(:,1)-180)<1e-3);
%phi270 = find(abs(PHI(:,1)-270)<1e-3);

% Co-polar and cx-polar phi=0 (H-plane)
Eco00 = [fliplr(Eco(phi180, 2:end)) Eco(phi000, :)];
Ecx00_t = [fliplr(Ecx(phi180, 2:end)) Ecx(phi000, :)];

% Co-polar and cx-polar phi=90 (E-plane)
Eco90 = [fliplr(Eco(phi270, 2:end)) Eco(phi090, :)];
Ecx90_t = [fliplr(Ecx(phi270, 2:end)) Ecx(phi090, :)];

% Directivity
%D = Directivity(Ethe, Ephi, vthe, vphi);
%D0 = 10*log10(max(max(D)));
u = (1/sqrt(2))*((sqrt(lambda*rho_2))/a1 + (a1/sqrt(lambda*rho_2)));
v = (1/sqrt(2))*((sqrt(lambda*rho_2))/a1 - (a1/sqrt(lambda*rho_2)));
Cu = fresnelc(u);
Cv = fresnelc(v);
Su = fresnels(u);
Sv = fresnels(v);
aux = b1/(sqrt(2*lambda*rho_1));
C2 = fresnelc(aux)^2;
S2 = fresnels(aux)^2;
Dp = (8*pi*rho_1*rho_2)/(a1*b1)*((Cu-Cv)^2 + (Su - Sv)^2)*(C2 + S2);
D0 = 10*log10(max(max(Dp)));

figure(); 
hold on; 
plot(vthe, D0+20*log10(abs(Eco00)), LineWidth=2)
plot(vthe, D0+20*log10(abs(Eco90)), LineWidth=2)
%plot(vthe, D0+20*log10(abs(Ecx00_t)), LineWidth=2)
%plot(vthe, D0+20*log10(abs(Ecx90_t)), LineWidth=2)
grid on; 
xlim([-90 90])
ylabel('Directivity [dB]', 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('\theta [º]', 'FontName', 'Times New Roman', 'FontSize', 14)
legend('Copolar H-Plane, \phi = 0º', 'Copolar E-Plane, \phi = 90º',...
    'FontName', 'Times New Roman', 'FontSize', 14)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

CbandTheoretycalCoHPlane= D0+20*log10(abs(Eco00));
save('Cband_theoretycal_Co_HPlane.mat', 'CbandTheoretycalCoHPlane');
CbandTheoretycalCoEPlane= D0+20*log10(abs(Eco90));
save('Cband_theoretycal_Co_EPlane.mat', 'CbandTheoretycalCoEPlane');

