clear,clc
%Set angle and distance inputs
phi = input('Enter the Latitude(Deg): ');
phi = phi * pi/180;
lambda = input('Enter the Local Sidereal Time(Deg): ');
lambda = lambda * pi / 180;
rho = input('Enter the Slant Range(km): ');
beta = input('Enter the Azimuth angle(Deg): ');
beta = beta * pi / 180; 
sigma = input('Enter the Elevation Angle(Deg): ');
sigma = sigma * pi / 180;
rho_dot = input('Enter the Range Rate(km/s): ');
beta_dot = input('Enter the angular rate of the Azimuth Angle(Deg/s): ');
beta_dot = beta_dot * pi /180 ;
sigma_dot = input('Enter the angular rate of the Elevation Angle(Deg/s): ');
sigma_dot = sigma_dot * pi / 180; 
TOF = input('Enter the Time of Flight for Case 1(min): ');
TOF = TOF * 60;
TOF2 = input('Enter the Time of Flight for Case 2(min): ');
TOF2 = TOF2 * 60;

%Run function OrbitComp with and without J2 effects for Case 1 and 2
J2 = false;
[r0, v0, OE0,rf,vf,OEf] = OrbitComp(phi, lambda, rho,beta, sigma, rho_dot, beta_dot, sigma_dot, TOF,J2);
[r02, v02, OE02,rf2,vf2,OEf2] = OrbitComp(phi, lambda, rho,beta, sigma, rho_dot, beta_dot, sigma_dot, TOF2,J2);

J2 = true;
[r0_J2, v0_J2, OE0_J2,rf_J2,vf_J2,OEf_J2] = OrbitComp(phi, lambda, rho,beta, sigma, rho_dot, beta_dot, sigma_dot, TOF,J2);
[r0_J2_2, v0_J2_2, OE0_J2_2,rf_J2_2,vf_J2_2,OEf_J2_2] = OrbitComp(phi, lambda, rho,beta, sigma, rho_dot, beta_dot, sigma_dot, TOF2,J2);


y0 = [r0; v0];

%Integration time span calcs
a = OE0(1,1);
u = 398600.44;
T = 2 * pi * sqrt(a^3/u);
tspan = [0 2*T];

%Integrate using ode45
tolerance = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t, y] = ode45(@(t, y) [y(4:6); -u / norm(y(1:3))^3 * y(1:3)], [0 2*T], [r0; v0],tolerance);
[tf, yf] = ode45(@(t, y) [y(4:6); -u / norm(y(1:3))^3 * y(1:3)], [0 TOF], [r0; v0],tolerance);
rf_ode45 = yf(end,1:3);
vf_ode45 = yf(end,4:6);

%Plot the orbit
plot3(y(:,1), y(:,2), y(:,3),'r','LineWidth',2);
hold on 
Re = 6378.1366;
[x, y, z] = sphere(100);
surf(Re*x,Re*y,Re*z,'FaceColor','blue','FaceAlpha',0.3)%Creates earth reference in plot;
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
grid on;
axis equal; 

%Display Results
disp('================== INITIAL POSITION AND VELOCITY ==================');
fprintf('Initial Position r0 [km]:     [%.10f  %.10f  %.10f]\n', r0);
fprintf('Initial Velocity v0 [km/s]:   [%.10f  %.10f  %.10f]\n', v0);

disp('================== FINAL STATE (OEtoOE) Case 1 ==================');
fprintf('Final Position rf [km]:       [%.10f  %.10f  %.10f]\n', rf);
fprintf('Final Velocity vf [km/s]:     [%.10f  %.10f  %.10f]\n', vf);

disp('================== FINAL STATE (OEtoOE) Case 2 ==================');
fprintf('Final Position rf [km]:       [%.10f  %.10f  %.10f]\n', rf2);
fprintf('Final Velocity vf [km/s]:     [%.10f  %.10f  %.10f]\n', vf2);

disp('================== FINAL STATE (ODE45) Case 1 ==================');
fprintf('Final Position rf [km]:       [%.10f  %.10f  %.10f]\n', rf_ode45);
fprintf('Final Velocity vf [km/s]:     [%.10f  %.10f  %.10f]\n', vf_ode45);

disp('================== INITIAL ORBITAL ELEMENTS ==================');
fprintf('a (km):        %.10f\n', OE0(1));
fprintf('e:             %.10f\n', OE0(2));
fprintf('i (deg):       %.10f\n', OE0(3));
fprintf('Omega (deg):  %.10f\n', OE0(4));
fprintf('omega (deg):   %.10f\n', OE0(5));
fprintf('f (deg):       %.10f\n', OE0(6));

disp('================== FINAL ORBITAL ELEMENTS Case 1 ==================');
fprintf('a (km):        %.10f\n', OEf(1));
fprintf('e:             %.10f\n', OEf(2));
fprintf('i (deg):       %.10f\n', OEf(3));
fprintf('Omega (deg):  %.10f\n', OEf(4));
fprintf('omega (deg):   %.10f\n', OEf(5));
fprintf('f (deg):       %.10f\n', OEf(6));

disp('================== FINAL ORBITAL ELEMENTS Case 2 ==================');
fprintf('a (km):        %.10f\n', OEf2(1));
fprintf('e:             %.10f\n', OEf2(2));
fprintf('i (deg):       %.10f\n', OEf2(3));
fprintf('Omega (deg):  %.10f\n', OEf2(4));
fprintf('omega (deg):   %.10f\n', OEf2(5));
fprintf('f (deg):       %.10f\n', OEf2(6));

disp('================== FINAL ORBITAL ELEMENTS W/ J2 Case 1 ==================');
fprintf('a (km):        %.10f\n', OEf_J2(1));
fprintf('e:             %.10f\n', OEf_J2(2));
fprintf('i (deg):       %.10f\n', OEf_J2(3));
fprintf('Omega (deg):  %.10f\n', OEf_J2(4));
fprintf('omega (deg):   %.10f\n', OEf_J2(5));
fprintf('f (deg):       %.10f\n', OEf_J2(6));

disp('================== FINAL ORBITAL ELEMENTS W/ J2 Case 2 ==================');
fprintf('a (km):        %.10f\n', OEf_J2_2(1));
fprintf('e:             %.10f\n', OEf_J2_2(2));
fprintf('i (deg):       %.10f\n', OEf_J2_2(3));
fprintf('Omega (deg):  %.10f\n', OEf_J2_2(4));
fprintf('omega (deg):   %.10f\n', OEf_J2_2(5));
fprintf('f (deg):       %.10f\n', OEf_J2_2(6));



