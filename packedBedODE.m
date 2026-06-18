function dydt = packedBedODE(~, y, n, dz, u, D_ax, eps, rho_s, k_a, k_d, q_max, C_in)
    % Split variables
    C = max(0, y(1:n));         % Gas concentration
    q = max(0, min(y(n+1:end), q_max));  % Adsorbed phase
    
    % Initialize derivatives
    dCdt = zeros(n,1);
    dqdt = zeros(n,1);
    
    % Boundary conditions
    C_left = C_in;
    C_right = C(end);  % Open boundary
    
    for i = 1:n
        % Upwind differencing for stability
        if i == 1
            dCdz = (C(1) - C_left)/dz;
            d2Cdz2 = (C(2) - 2*C(1) + C_left)/dz^2;
        elseif i == n
            dCdz = (C(n) - C(n-1))/dz;
            d2Cdz2 = (C_right - 2*C(n) + C(n-1))/dz^2;
        else
            dCdz = (C(i) - C(i-1))/dz;
            d2Cdz2 = (C(i+1) - 2*C(i) + C(i-1))/dz^2;
        end
        
        % Adsorption kinetics with limits
        q_avail = max(0, q_max - q(i));
        dqdt(i) = min(0.1, max(-0.1, k_a*C(i)*q_avail - k_d*q(i)));
        
        % Mass balance
        dCdt(i) = (1/eps)*(D_ax*d2Cdz2 - u*dCdz) - (1-eps)/eps*rho_s*dqdt(i);
    end
    
    dydt = [dCdt; dqdt];
end