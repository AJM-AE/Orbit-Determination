function f2 = OEtoOE(a,e,delta_t,M1)
u = 398600;
n = sqrt(u/a^3);
M2 = M1 + n*delta_t;

E2 = M2;

converge = false;
tol = 1e-5;
while converge == false
    f = E2 - e * sin(E2) - M2;  % Kepler's equation
    f_prime = 1 - e * cos(E2); % Derivative of f(E)

    E_new = E2 - f / f_prime;

    % Check for convergence
    if abs(E_new - E2) < tol
        converge = true;
    end
    
    E2 = E_new;
end

f2 = 2*atan(sqrt((1+e)/(1-e))*tan(E2/2));
