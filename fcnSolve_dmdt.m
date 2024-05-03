
function [t,rho, Tair, P, Trock,r,Nc] = fcnSolve_dmdt(NbCycles,para,SolverOptions,flagRealGas)
%
%   This function will solve the differential coupled system of equations for
%   the temperature, pressure and density of hydrogen gaz injected in
%   a tank (outside) or a cavity (in the ground) 
%   
%   Heat exchange coefficient can be imposed or will be estimated from Nusslet parametrs
%
%   A function for dm/dt is imposed 
%   during charge and discharge. (Real gaz model is used for hydrogen with
%   a compressibility empirical model) 
%
%   THE INPUTS:
%       
%       1- NbCycles 
%               The code solve the differential system of equations for "NbCycles" identical  
%               CAES plant cycle (NbCycles HAS TO BE AN INTEGER). A CAES plant cycle
%               starts at t=0 and finiches at t = 1 (t is a normalized time)
%
%       2- para
%               The structure "para" is used to pass all the normalized parameters 
%               and the interpolation functions for the hydrogen properties  
%               This strucutre is set with the function ParametersForSolver(p)
%               called from the main file.
%
%       3- SolverOptions
%               A structure for the options of the numerical integration
%               This strucutre is set in the main file.
%
%       4- flgagRealGas
%               1- is the real hydrogen gaz with compressibility given by
%               Lenon equation
%               0- is a simplifed gaz (constant compressibility or
%               compressibility equal to one) depending on the choice made
%               in the main file
%
%   The coupled system of equations is store in variable Y:
%   Y[1] -> Temperature in the rock at r = 1 (j=1) (at cavern wall)
%   Y[2] -> Temperature in the rock at r = r (j=2) (eq 24 of Kushnir; j starts at 1 in this code!)
%     :
%   Y[N+1] -> Temperature in the rock at r = Rp
%   Y[N+2] -> air temperature (Normalized)
%   Y[N+3] -> density (Normalized)
%
%   THE FUNCTION RETURNS:
%
%   t   : the time where the solution was evaluated (normalized time)
%   rho : the normalized air density at different time
%   Tair: the normalized air temperature
%   Trock: the normalized rock temperature at r (N+1 values of r: r=1..r=Rp)
%   ml: normalized leakage rate (set to 0 so far...)
%   Nc: index to indicate the end of cycle 
%           Nc(1) -> t(Nc(1)) = 1
%           Nc(2) -> t(Nc(2)) = 2
%             :
%           Nc(NbCycles) -> t(Nc(NbCycles)) = NbCycles 
% 
% ======================================================================
% Copyright (c) June 2023, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Numerical options for the grid size  
    N = SolverOptions.N;    
    beta = SolverOptions.beta;
    flagRp = floor(SolverOptions.Rp);
    Rp = SolverOptions.Rp;
    eps = SolverOptions.eps;
    
    if (SolverOptions.solver)
        solver = @ode23s;
    else
        solver = @ode45;
    end
    options = odeset('RelTol',SolverOptions.RelTol,'AbsTol',SolverOptions.AbsTol);
    NtimeStep = SolverOptions.NtimeStep;
         
    % initial condition
    Yinitial = ones(1,N+1)*para.TRw0; % Y(1:N+1): Rock temperature at N+1 values of r 
    Yinitial = [Yinitial 1 1];          % Y(N+2):   Air temp
                                        % Y(N+3):   rho
                                 
    eta = linspace(0,1,N+1);    % eta parameter from 0 to 1 (N+1 values)
    eta_p = zeros(1,N+1);       % first derivative of eta with respect to r
    eta_pp = zeros(1,N+1);      % second derivative of eta with respect to r
    
    kR = zeros(1,N+1); 
    F0 = zeros(1,N+1);
    Nc = zeros(1,NbCycles);  
    
    if (para.underground == 0)
        Rp = para.Rw(end);
        r = fcnrj(eta,beta,Rp);
        Rchange = para.Rw(1);
        cnt = 1;
        for j = 1:N+1
            eta_p(j) = fcneta_p(r(j),beta,Rp);    % calculate first derivative of eta
            eta_pp(j) = fcneta_pp(r(j),beta,Rp);  %  calculate second derivative of eta
            if (r(j)>Rchange)
                kR(j) = para.kR(cnt+1);
                F0(j) = para.F0(cnt+1);
                Rchange = para.Rw(cnt+1);
                cnt = cnt+1;
            else
                kR(j)= para.kR(cnt);
                F0(j) = para.F0(cnt);
            end
        end
    end   
    
    % main loop
    maxFokR = max(para.F0.*para.kR);
    for n = 1:NbCycles
        if (para.underground)
            % Evaluate Rp from eq A.7 of Kushnir
            if not(flagRp)
                xi_n = sqrt(2*pi/(maxFokR*n));
                Rp = 1+sqrt(2)/xi_n*log(1/eps);
            end
            r(n,:) = fcnrj(eta,beta,Rp);
            Rchange = para.Rw(1);
            cnt = 1;
            for j = 1:N+1
                eta_p(j) = fcneta_p(r(n,j),beta,Rp);    % calculate first derivative of eta
                eta_pp(j) = fcneta_pp(r(n,j),beta,Rp);  %  calculate second derivative of eta
                if (r(n,j)>Rchange)
                    kR(j) = para.kR(cnt+1);
                    F0(j) = para.F0(cnt+1);
                    Rchange = para.Rw(cnt+1);
                    cnt = cnt+1;
                else
                    kR(j)= para.kR(cnt);
                    F0(j) = para.F0(cnt);
                end
            end
        end
        tspan = [0 1];
        switch (flagRealGas)
            case 0 % Simple or ideal gas model 
                if (NtimeStep) % We read the solution at NtimeStep + 1 (t=0)
                    sol = solver(@(t,Y) odefcn_dmdt_Simple(t,Y,r(n,:),kR,F0,eta_p,eta_pp,para),tspan,Yinitial,options);
                    tspan = linspace(0,1,NtimeStep+1);
                    Yspan = (deval(sol,tspan)).';
                else % we take only the time steps choosen by the numerical integrator
                    [tspan,Yspan] = solver(@(t,Y) odefcn_dmdt_Simple(t,Y,r(n,:),kR,F0,eta_p,eta_pp,para),tspan,Yinitial,options);
                end
            case 1 % Real gas model 
                 if (NtimeStep) % We read the solution at NtimeStep + 1 (t=0)
                    sol = solver(@(t,Y) odefcn_dmdt_Real(t,Y,r(n,:),kR,F0,eta_p,eta_pp,para),tspan,Yinitial,options);
                    tspan = linspace(0,1,NtimeStep+1);
                    Yspan = (deval(sol,tspan)).';
                else % we take only the time steps choosen by the numerical integrator
                    [tspan,Yspan] = solver(@(t,Y) odefcn_dmdt_Real(t,Y,r(n,:),kR,F0,eta_p,eta_pp,para),tspan,Yinitial,options);
                end
        end
        Yinitial = Yspan(end,:); % set the inttial condition for next time interval
            % store the calculated values
            % the number of calculated values will change from cycle to
            % cycle
         if (n==1)
             t = (n-1)+tspan;
             rho = Yspan(:,N+3);
             Z = para.Z_Trho(Yspan(:,N+2)*para.T0,Yspan(:,N+3)*para.rho0);
             P =  Yspan(:,N+3).*Yspan(:,N+2).*Z/para.Z0;
             Tair = Yspan(:,N+2);
             Trock = Yspan(:,1:N+1);
         else
             t = [t;(n-1)+tspan]; %#ok<*AGROW>
             rho = [rho;Yspan(:,N+3);];
             Z = para.Z_Trho(Yspan(:,N+2)*para.T0,Yspan(:,N+3)*para.rho0);
             P = [P;Yspan(:,N+3).*Yspan(:,N+2).*Z/para.Z0];
             Tair = [Tair;Yspan(:,N+2)];
             Trock = [Trock;Yspan(:,1:N+1)];
         end
         Nc(n) = Nc(n) + length(tspan);
    end
end

function dYdt = odefcn_dmdt_Real(t,Y,r,kR,F0,eta_p,eta_pp,para)

    % In this function, the density (Y(N+3)) and temperature (Y(N+2))  
    % are numerically intergrated. The pressure Y(N+4) is calculated from
    % the density, the temperature and de compressibility Z. Hence, the
    % pressure is not coupled with the differential equations system.
   
    % In this function gamma is only a flag (0 or 1) 
    
    % Charge: Fie >0 AND gamma = 1
    % Discharge Fie <0 AND gamma = 0
    % Storage Fie = 0
    
    if (para.Fie(t)>0)
        gamma = 1;
    else
        gamma = 0;
    end
    
    N = length(Y)-3;
    dYdt = zeros(N+3,1);
    
    T = Y(N+2)*para.T0; % temperature in oK
    rho = Y(N+3)*para.rho0; % density in kg/m^3
    Z = para.Z_Trho(T,rho);
    P = Y(N+3)*Y(N+2)*Z/para.Z0*para.P0;
    Ti = para.Ti*para.T0;
 
    dZ = para.dZdT_rho(T,P)*para.T0;
    Hi = para.H(Ti,P);
    Hc = para.H(T,P);
    cv = para.Cv(T,P);
    
    beta = 1/Y(N+2) + para.T0/Z*para.dZdT_P(T,P);
    
    if (para.Nusselt) % try to estimate heat exchange coefficient from Nusselt parameter
        
        mdot = para.mr/para.tp*para.rho0*para.V;
        
        %density of injected gaz;
        rho_i = Z*rho*T/para.Z_TP(Ti,P)/Ti;
        % velocity of injected gaz
        v = mdot./rho_i./pi/(para.Dinlet/2)^2;
        if (v > 50)
            v = 50;
        end
        
        mu_inlet = para.mu(Ti,P);
        mu = para.mu(T,P);
        kH2 = para.kH2(T,P);
        cp = para.Cp(T,P);
        Beta = beta/para.T0;     % beta is beta* in my notes
        
        Re_in = 4*mdot/mu_inlet/pi/para.Dinlet;
        Ra = 9.8*Beta*cp*rho^2*para.L^3*abs(T-Y(1)*para.T0)/mu/kH2;
        
        a = 0.323;
        Nu_for = a*Re_in^0.67;
        Nu_nat = 0.104*Ra^0.352;
        hc = kH2*(Nu_for + Nu_nat)/para.L;
        Bi_int = hc*para.Bi_int;
        qr = hc*para.qr;
    else
        Bi_int = para.Bi_int;
        qr = para.qr;
        mdot = para.mr/para.tp*para.rho0*para.V;
        %density of injected gaz;
        rho_i = Z*rho*T/para.Z_TP(Ti,P)/Ti;
        % velocity of injected gaz
        v = mdot./rho_i./pi/(para.Dinlet/2)^2;
        if (v > 50)
            v = 50;
        end
    end
   
    ml = 0;
    
    kRp = (kR(2)+kR(1))/2;
    dYdt(1) = F0(1)*kRp*(2*N^2*eta_p(1)^2*(Y(2)-Y(1))) + F0(1)*kR(1)*(Bi_int*(eta_pp(1)/eta_p(1) + 1/r(1) - 2*N*eta_p(1))*(Y(1)-Y(N+2)));
      
    for j=2:N
       kRp = (kR(j+1)+kR(j))/2;
       kRm = (kR(j-1)+kR(j))/2;
       dYdt(j) = F0(j)*(N^2*eta_p(j)^2*(kRp*(Y(j+1)-Y(j))+kRm*(Y(j-1)-Y(j)))+ N*kR(j)*(eta_pp(j)+eta_p(j)/r(j))*(Y(j+1)-Y(j-1))/2);
    end
    kRm = (kR(N+1)+kR(N))/2;
    if (para.underground)
        dYdt(N+1) = F0(N+1)*kRm*(2*N^2*eta_p(N+1)^2*(Y(N)-Y(N+1)));
    else
        dYdt(N+1) = F0(N+1)*kRm*(2*N^2*eta_p(N+1)^2*(Y(N)-Y(N+1))) + F0(N+1)*kR(N+1)*(para.Bi_ext*(eta_pp(N+1)/eta_p(N+1) + 1/r(N+1) + 2*N*eta_p(N+1))*(para.Tair-Y(N+1)));
    end
    
    dYdt(N+2) = para.mr/Y(N+3)/cv*(para.Fie(t)*gamma*(Hi-Hc)/para.T0 + (para.Fie(t)-ml)*(Z*para.R*Y(N+2) + para.R*(Y(N+2))^2*dZ + para.potential/para.T0 + v^2/para.T0) + qr*(Y(1)-Y(N+2)));    
    dYdt(N+3) = (para.Fie(t))*para.mr;
end


function dYdt = odefcn_dmdt_Simple(t,Y,r,kR,F0,eta_p,eta_pp,para)
   
    % Charge: Fie >0
    % Discharge Fie <0 AND gamma = 0
    % Storage Fie = 0
    
    if (para.Fie(t)<0)
        gamma = 0;
    else
        gamma = para.gamma;
    end
   
    N = length(Y)-3;
    dYdt = zeros(N+3,1);
    
    T = Y(N+2)*para.T0; % temperature in oK
    rho = Y(N+3)*para.rho0; % density in kg/m^3
    Z = para.Z_Trho(T,rho);
    P = Y(N+3)*Y(N+2)*Z/para.Z0;
    Ti = para.Ti*para.T0;

    if (para.ZT0==0)
        beta = 1/Y(N+2);    %Ideal gaz
    else
        beta = 1/Y(N+2) + para.T0/para.Z0*para.dZdT_P(para.T0,para.P0); % Simpliifed model 
    end  
    
    ml = 0;
    
    if (para.Nusselt) % try to estimate heat exchange coefficient from Nusselt parameter
        
        mdot = para.mr/para.tp*para.rho0*para.V;
        
        %density of injected gaz;
        rho_i = Z*rho*T/para.Z_TP(Ti,P)/Ti;
        % velocity of injected gaz
        v = mdot./rho_i./pi/(para.Dinlet/2)^2;
        if (v > 50)
            v = 50;
        end
        
        mu_inlet = para.mu(Ti,P);
        mu = para.mu(T,P);
        kH2 = para.kH2(T,P);
        cp = para.Cp(T,P);
        Beta = beta/para.T0;     % beta is beta* in my notes
        
        Re_in = 4*mdot/mu_inlet/pi/para.Dinlet;
        Ra = 9.8*Beta*cp*rho^2*para.L^3*abs(T-Y(1)*para.T0)/mu/kH2;
        
        a = 0.323;
        Nu_for = a*Re_in^0.67;
        Nu_nat = 0.104*Ra^0.352;
        hc = kH2*(Nu_for + Nu_nat)/para.L;
        Bi_int = hc*para.Bi_int;
        qr = hc*para.qr;
    else
        Bi_int = para.Bi_int;
        qr = para.qr;
        %{
        mdot = para.mr/para.tp*para.rho0*para.V;
        %density of injected gaz;
        rho_i = Z*rho*T/para.Z_TP(Ti,P)/Ti;
        % velocity of injected gaz
        v = mdot./rho_i./pi/(para.Dinlet/2)^2;
        %}
      %  v = 0; % to compare with analytical solution
    end
    
    kRp = (kR(2)+kR(1))/2;
    dYdt(1) = F0(1)*kRp*(2*N^2*eta_p(1)^2*(Y(2)-Y(1))) + F0(1)*kR(1)*(Bi_int*(eta_pp(1)/eta_p(1) + 1/r(1) - 2*N*eta_p(1))*(Y(1)-Y(N+2)));
   
    for j=2:N
        kRp = (kR(j+1)+kR(j))/2;
        kRm = (kR(j-1)+kR(j))/2;
        dYdt(j) = F0(j)*(N^2*eta_p(j)^2*(kRp*(Y(j+1)-Y(j))+kRm*(Y(j-1)-Y(j)))+ N*kR(j)*(eta_pp(j)+eta_p(j)/r(j))*(Y(j+1)-Y(j-1))/2);
    end
    kRm = (kR(N+1)+kR(N))/2;
    if (para.underground)
        dYdt(N+1) = F0(N+1)*kRm*(2*N^2*eta_p(N+1)^2*(Y(N)-Y(N+1)));
    else
        dYdt(N+1) = F0(N+1)*kRm*(2*N^2*eta_p(N+1)^2*(Y(N)-Y(N+1))) + F0(N+1)*kR(N+1)*(para.Bi_ext*(eta_pp(N+1)/eta_p(N+1) + 1/r(N+1) + 2*N*eta_p(N+1))*(para.Tair-Y(N+1)));
    end
    
    dYdt(N+2) = para.mr/Y(N+3)*(para.Fie(t)*gamma*(para.Ti-Y(N+2))+(para.Fie(t)-ml)*(para.Rs*Y(N+2)+para.Us*(Y(N+2))^2) + qr/para.cv0*(Y(1)-Y(N+2)));
    dYdt(N+3) = (para.Fie(t)-ml)*para.mr;
end



function rj = fcnrj(eta,beta,Rp)
    % the values of r for the uniform interval of eta (0<eta<1) 
    rj = 1+(Rp-1)*(beta+1-(beta-1)*((beta+1)/(beta-1)).^(1-eta))./(1+((beta+1)/(beta-1)).^(1-eta));
end

function eta_p = fcneta_p(r,beta,Rp)
    % first derivative of eta with respect to r
  eta_p = -(-1/(Rp-1)/(beta-1+(r-1)/(Rp-1))-(beta+1-(r-1)/(Rp-1))/(beta-1+(r-1)/(Rp-1))^2/(Rp-1))/(beta+1-(r-1)/(Rp-1))*(beta-1+(r-1)/(Rp-1))/log(((beta+1)/(beta-1)));
end

function eta_pp = fcneta_pp(r,beta,Rp)
    % second derivative of eta with respect to r
   eta_pp = -(2 / (Rp - 1) ^ 2 / (beta - 1 + (r - 1) / (Rp - 1)) ^ 2 + 2 * (beta + 1 - (r - 1) / (Rp - 1)) / (beta - 1 + (r - 1) / (Rp - 1)) ^ 3 / (Rp - 1) ^ 2) / (beta + 1 - (r - 1) / (Rp - 1)) * (beta - 1 + (r - 1) / (Rp - 1)) / log(((beta + 1) / (beta - 1))) - (-1 / (Rp - 1) / (beta - 1 + (r - 1) / (Rp - 1)) - (beta + 1 - (r - 1) / (Rp - 1)) / (beta - 1 + (r - 1) / (Rp - 1)) ^ 2 / (Rp - 1)) / ((beta + 1 - (r - 1) / (Rp - 1)) ^ 2) * (beta - 1 + (r - 1) / (Rp - 1)) / log(((beta + 1) / (beta - 1))) / (Rp - 1) - (-1 / (Rp - 1) / (beta - 1 + (r - 1) / (Rp - 1)) - (beta + 1 - (r - 1) / (Rp - 1)) / (beta - 1 + (r - 1) / (Rp - 1)) ^ 2 / (Rp - 1)) / (beta + 1 - (r - 1) / (Rp - 1)) / (Rp - 1) / log(((beta + 1) / (beta - 1)));
end



