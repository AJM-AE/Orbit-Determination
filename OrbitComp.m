function [r0, v0, OE0,rf,vf,OEf] = OrbitComp(phi, lambda, rho, beta, sigma, rho_dot, beta_dot, sigma_dot, TOF,J2)

%Sets vectors for angular velocity of earth and the site location in SEZ
W_E = [0 ; 0 ; 7.2921 * 10^-5]; 
SEZ_rsite = [0;0;6378];

%Calculates the position and velocity vectors in SEZ
rho_SEZ = [-rho*cos(sigma)*cos(beta); rho*cos(sigma)*sin(beta);rho*sin(sigma)];

rho_dot_SEZ = [-rho_dot*cos(sigma)*cos(beta) + rho*sigma_dot*sin(sigma)*cos(beta) + rho*beta_dot*cos(sigma)*sin(beta); 
    rho_dot*cos(sigma)*sin(beta) - rho*sigma_dot*sin(sigma)*sin(beta) + rho*beta_dot*cos(sigma)*cos(beta);
    rho_dot*sin(sigma) + rho*sigma_dot*cos(sigma)];

%Converts position and velocity vectors to ECI
RotT = [cos(lambda) -sin(lambda) 0 ; sin(lambda) cos(lambda) 0; 0 0 1] * [sin(phi) 0 cos(phi); 0 1 0; -cos(phi) 0 sin(phi)];

rho_ECI = RotT * rho_SEZ;

ECI_rsite = RotT * SEZ_rsite;

r0 = rho_ECI + ECI_rsite;

rho_dot_ECI = RotT * rho_dot_SEZ;

v0 = rho_dot_ECI + cross(W_E,rho_ECI) + cross(W_E,ECI_rsite);

I = [1;0;0];
J = [0;1;0];
K = [0;0;1];

%Calculates the size and shape orbital elements 

u = 398600.44;
h = cross(r0,v0);
h_unit = h/norm(h);
h_magnitude = norm(h);

r_unit = r0/norm(r0);
r_magnitude = norm(r0);

v_magnitude = norm(v0);

epsilon = (v_magnitude^2)/2 - u/r_magnitude;

a = -u/(2*epsilon);

e = cross(v0,h)/398600 - r_unit;
e_unit = e/norm(e);
e_magnitude = norm(e);

%Calculates orientation orbital elements
i = acos(dot(K,h_unit));

n_unit = cross(K, h_unit) / norm(cross(K, h_unit));

Omega = atan2(n_unit(2), n_unit(1));

quadrant_check1 = dot(e_unit, K);
if quadrant_check1 >= 0
    omega = acos(dot(n_unit, e_unit));
else
    omega = 2*pi - acos(dot(n_unit, e_unit));
end

%Calculates f at t1
quadrant_check2 = dot(r0,v0);
if quadrant_check2 >= 0
    f = acos(dot(e_unit, r_unit));
else
    f = 2*pi - acos(dot(e_unit, r_unit));
end

% Eccentric anomaly calculations
if e_magnitude < 1  % Elliptical orbit
    E1 = 2 * atan(sqrt((1 - e_magnitude) / (1 + e_magnitude)) * tan(f / 2));
    M1 = E1 - e_magnitude * sin(E1);
    n = sqrt(u / a^3);
    t1 = M1 / n;

elseif e_magnitude > 1  % Hyperbolic orbit 
    sinhF = sqrt(e_magnitude^2 - 1) * sin(f) / (1 + e_magnitude * cos(f));
    F = asinh(sinhF); 
    M1 = e_magnitude * sinhF - F;
    n = sqrt(u / (-a)^3);
    t1 = M1 / n;
end

%Converts OE angles to degrees
i_deg = i / pi * 180;
Omega_deg = Omega / pi * 180;
omega_deg = omega / pi * 180;
f_deg = f / pi * 180;

%Creates an OE array for t1
OE0 = [a;e_magnitude;i_deg;Omega_deg;omega_deg;f_deg];

%Calculates f at t2
f2 = OEtoOE(a,e_magnitude,TOF,M1);

%Adjusts Omega and omega with J2
if J2 == true
    n = sqrt(u / a^3);
    j2 = 1.08264 * 10^-3;
    Re = 6378.1366;
    Omega_dot = (-3*n*j2 / (2*(1-e_magnitude^2)^2))*(Re/a)^2 * cos(i);
    omega_dot = (3*n*j2 / (4*(1-e_magnitude^2)^2))*(Re/a)^2 * (5*cos(i)^2 -1);
    Omega = Omega + Omega_dot * TOF;
    omega = omega + omega_dot * TOF; 
    Omega_deg = Omega / pi * 180;
    omega_deg = omega / pi * 180;
end

%Converts f2 to degrees and creates an array for the final OEs
f2_deg = f2 / pi * 180;
OEf = [a;e_magnitude;i_deg;Omega_deg;omega_deg;f2_deg];

%Converts OEf to rf and vf

u = 398600;
p = a*(1-e_magnitude^2);
r2_P = [p/(1 + e_magnitude * cos(f2)) * cos(f2);p/(1 + e_magnitude*cos(f2))*sin(f2);0];
v2_P = [u/h_magnitude*-sin(f2); u/h_magnitude*(e_magnitude + cos(f2));0];

rot3_Omega = [cos(-Omega), sin(-Omega), 0; -sin(-Omega), cos(-Omega), 0; 0, 0, 1];
rot1_i = [1, 0, 0; 0, cos(-i), sin(-i); 0, -sin(-i), cos(-i)];
rot3_omega = [cos(-omega), sin(-omega), 0; -sin(-omega), cos(-omega), 0; 0, 0, 1];

rf = rot3_Omega * rot1_i * rot3_omega * r2_P;
vf = rot3_Omega * rot1_i * rot3_omega * v2_P;

