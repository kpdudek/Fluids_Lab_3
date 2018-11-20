function FluidsLab3
load('Pressure_Data.mat')
a = .01; %m
inletP = [.11132,1.1104,.00035556,.026633]; %Pa
u_oil = 10;
u_water = 8.90 * 10^-4;
U_oil_1 = .00051;
U_oil_2 = .0051;
U_water_1 = .005;
U_water_2 = .05;
p_water = 997;
p_oil = 971;

pressure_vs_theta_OIL(surface_pressure_SiliconOilRe001,inletP(1),a,u_oil,U_oil_1)
pressure_vs_theta_OIL(surface_pressure_SiliconOilRe2,inletP(2),a,u_oil,U_oil_2)

pressure_vs_theta_WATER(surface_pressure_WaterRe1,inletP(3),a,p_water,U_water_1)
pressure_vs_theta_WATER(surface_pressure_WaterRe2,inletP(4),a,p_water,U_water_2)

D = [.0011712,.011712,4.3567*10^-6,1.7319*10^-4];
[cdO1,ReO1,cdTO1] = calculate_Re_and_Cd(D(1),p_oil,u_oil,a,U_oil_1);
[cdO2,ReO2,cdTO2] = calculate_Re_and_Cd(D(2),p_oil,u_oil,a,U_oil_2);
[cdW1,ReW1,cdTW1] = calculate_Re_and_Cd(D(3),p_water,u_water,a,U_water_1);
[cdW2,ReW2,cdTW2] = calculate_Re_and_Cd(D(4),p_water,u_water,a,U_water_2);

Re = [ReO1,ReO2,ReW1,ReW2];
cd = [cdO1,cdO2,cdW1,cdW2];
cdT = [cdTO1,cdTO2,cdTW1,cdTW2];

plot_drags(Re,cd,cdT)

function pressure_vs_theta_OIL(pressure,inlet,a,u,U)
delta_p = pressure(:,2) - inlet;
thetas = pressure(:,1)./a;

theoretical = (3/2) * cos(thetas);
dimensionless_P = delta_p / ((u*U)/a);

figure
plot(thetas,dimensionless_P,'b.')
hold on
plot(thetas,theoretical)

function pressure_vs_theta_WATER(pressure,inlet,a,p,U)
delta_p = pressure(:,2) - inlet;
thetas = pressure(:,1)./a;

dimensionless_P = delta_p / (.5*p*U^2);
theoretical = ones(length(pressure),1)-(((9/4).*((sin(thetas)).^2)));

figure
plot(thetas,dimensionless_P,'b.')
hold on
plot(thetas,theoretical)

function [cd,Re,cd_theory] = calculate_Re_and_Cd(D,p,u,a,U)
cd = D / (.5 * p * U^2 * pi * a^2);
Re = (2 * p * U * a) / u;
cd_theory = 24/Re;

function plot_drags(re,cd,cdT)
figure
loglog(re(1),cd(1),'.',re(2),cd(2),'.',re(3),cd(3),'.',re(4),cd(4),'.','MarkerSize',20)
hold on 
loglog(re,cdT,'LineWidth',1.5)

















