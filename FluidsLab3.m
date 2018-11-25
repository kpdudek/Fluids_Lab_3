function FluidsLab3
close all
clear
clc

load('Pressure_Data.mat')

% Inlet pressure from COMSOL
inletP = [.11132,1.1104,.00035556,.026633]; %Pa

% Drag force from COMSOL
D = [.0011712,.011712,4.3567*10^-6,1.7319*10^-4]; % N

% Radius of the sphere (m)
a = .01;

% Densitites of the fluids (kg/m^3)
p_water = 997;
p_oil = 971;

% Dynamic viscosities of the fluids (kg*s/m^2)
u_oil = 10;
u_water = 8.90 * 10^-4;

% Velocitites of the flows (m/s)
U_oil_1 = .00051;
U_oil_2 = .0051;
U_water_1 = .005;
U_water_2 = .05;


%%% Calculate the dimensionless pressure and plot against the theta value
[O1_theta,O1_dimensionless_P,O1_theoretical] = pressure_vs_theta_OIL(surface_pressure_SiliconOilRe001,inletP(1),a,u_oil,U_oil_1);
[O2_theta,O2_dimensionless_P,O2_theoretical] = pressure_vs_theta_OIL(surface_pressure_SiliconOilRe2,inletP(2),a,u_oil,U_oil_2);

[W1_theta,W1_dimensionless_P,W1_theoretical] = pressure_vs_theta_WATER(surface_pressure_WaterRe1,inletP(3),a,p_water,U_water_1);
[W2_theta,W2_dimensionless_P,W2_theoretical] = pressure_vs_theta_WATER(surface_pressure_WaterRe2,inletP(4),a,p_water,U_water_2);

thetas = [O1_theta,O2_theta,W1_theta,W2_theta];
dimensionless_P = [O1_dimensionless_P,O2_dimensionless_P,W1_dimensionless_P,W2_dimensionless_P];
theoretical_P = [O1_theoretical,O2_theoretical,W1_theoretical,W2_theoretical];

plot_pressure_vs_theta(dimensionless_P,theoretical_P,thetas)


%%% Calculste the Reynolds number of the flow and the coefficients of drag
[cdO1,ReO1,cdTO1] = calculate_Re_and_Cd(D(1),p_oil,u_oil,a,U_oil_1,'Oil 1');
[cdO2,ReO2,cdTO2] = calculate_Re_and_Cd(D(2),p_oil,u_oil,a,U_oil_2,'Oil 2');
[cdW1,ReW1,cdTW1] = calculate_Re_and_Cd(D(3),p_water,u_water,a,U_water_1,'Water 1');
[cdW2,ReW2,cdTW2] = calculate_Re_and_Cd(D(4),p_water,u_water,a,U_water_2,'Water 2');

% Store values in vectors
Re = [ReO1,ReO2,ReW1,ReW2];
cd = [cdO1,cdO2,cdW1,cdW2];
cdT = [cdTO1,cdTO2,cdTW1,cdTW2];

%Plot reynolds number versus velocity
plot_drags(Re,cd,cdT)

function [thetas,dimensionless_P,theoretical] =  pressure_vs_theta_OIL(pressure,inlet,a,u,U)
delta_p = pressure(:,2) - inlet;
thetas = pressure(:,1)./a;

theoretical = (3/2) * cos(thetas);
dimensionless_P = delta_p / ((u*U)/a);

function [thetas,dimensionless_P,theoretical] = pressure_vs_theta_WATER(pressure,inlet,a,p,U)
delta_p = pressure(:,2) - inlet;
thetas = pressure(:,1)./a;

dimensionless_P = delta_p / (.5*p*U^2);
theoretical = ones(length(pressure),1)-(((9/4).*((sin(thetas)).^2)));

function plot_pressure_vs_theta(dimensionless_P,theoretical_P,thetas)
%%% Silicon oil
figure('Name','Silicon Oil')

plot(thetas(:,1),dimensionless_P(:,1),'b.',thetas(:,2),dimensionless_P(:,2),'g.')
hold on
plot(thetas(:,1),theoretical_P(:,1),'r')

legend('Re = 9.904e^{-4}','Re = 9.904e^{-3}','Theoretical Value','Location','northeast')
title('$\displaystyle\frac{\Delta P}{\mu U \frac{1}{a}}$ vs $\theta$ for Silicon Oil','interpreter','latex')
ylabel('$\displaystyle\frac{\Delta P}{\mu U \frac{1}{a}}$','interpreter','latex')
xlabel('$\theta$','interpreter','latex')

%%% Water
figure('Name','Water')
plot(thetas(:,3),dimensionless_P(:,3),'b.',thetas(:,4),dimensionless_P(:,4),'g.')
hold on
plot(thetas(:,3),theoretical_P(:,3),'r')

legend('Re = 1.120e^{2}','Re = 1.120e^{3}','Theoretical Value','Location','northeast')
title('$\displaystyle\frac{\Delta P}{\mu U \frac{1}{a}}$ vs $\theta$ for Water','interpreter','latex')
ylabel('$\displaystyle\frac{\Delta P}{\mu U \frac{1}{a}}$','interpreter','latex')
xlabel('$\theta$','interpreter','latex')

function [cd,Re,cd_theory] = calculate_Re_and_Cd(D,p,u,a,U,NAME)
cd = D / (.5 * p * U^2 * pi * a^2);
Re = (2 * p * U * a) / u;
cd_theory = 24/Re;
fprintf('The Reynolds number for %s is: %.3e\n',NAME,Re)

function plot_drags(re,cd,cdT)
figure

loglog(re(1),cd(1),'.',re(2),cd(2),'.',re(3),cd(3),'.',re(4),cd(4),'.','MarkerSize',20)
hold on 
loglog(re,cdT,'LineWidth',1.5)

legend('Oil: Re = 9.904e^{-4}','Oil: Re = 9.904e^{-3}','Water: Re = 1.120e^{2}','Water: Re = 1.120e^{3}','Theoretical C_d','Location','northeast')
ylabel('C_d')
xlabel('Re')
title('C_d vs Re for Silicon Oil and Water')
grid on


















