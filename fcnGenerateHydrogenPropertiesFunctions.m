
function fcnGenerateHydrogenPropertiesFunctions()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file creates numerical functions (Interpolation) for the thermophysical
% properties of the Hydrogen gaz (and other properties like viscosity and
% thermal conductivity)
%
% The interpolation functions are created from the emperical Lemmon equation 
% that gives the compressibility factor Z as a function of T (Temperature) 
% and P (pressure) =>  Z(T,P)
% 
%  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4652867/
%   
% From the Lemmon Z(T,P) it is possible to calculate analytically
% dZ/dT and d^2Z/dT^2 at constant density. 
%
% Since the Lemmon equation is a sommation and using the fact the the functions 
% Z(T,P); dZdT and d2Z/dT2 are smooths, it is more efficient to generate a
% scatterInterpolant from multiple data generated from these functions. 
% The scatterInterpolant functions are saved in the file
% "ScatterInterpFunctions"
%
% Enthalpy, specific heats, thermal expansion and compressibility can be calculated from
% the compressibility and the gas law. 
%
% For conducitvity interpolation function, we used Assael and al. paper ( https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=907928 (Accessed September 29, 2023) ) 
%
% For the viscosity interpolation function, we used Muzny and al. paper (Journal of Chemical Engeneering data 2013, 58 pp 969-979
% dx.doi.org/10.1021/je301273j )
%
%   10 functions are saved
% 
%   1) Z_TP Returns Z at a given Temperature (K) and pressure (Pa)
%       ex: Z_TP(300,1e6) returns the compressibility at T =300 K and P =
%       1e6 Pa
%
%   2) Z_Trho : Returns Z at a given Temperature (K) and density (kg/m3)
%       ex: Z_Trho(300,10) returns the compressibility at T =300 K and rho = 10 kg/m^3
%
%   3) dZdT_rho : Returns dZdT (at constant rho) at a given Temperature (K) and Pressure (Pa)
%       ex: dZdT_rho(300,1e6) returns the derivative of the compressibility according to T (at constant rho) 
%                           for T =300 K and P = 1e6 Pa
%
%   4) dZdP_T:  Returns dZdP_T at a given Temperature (K) and Pressure (Pa)
%               at constant T
%
%   5) dZdT_P:   Returns dZdT_P (at constant Pressure) at a given Temperature (K) and Pressure (Pa)
%       ex: dZdTconstP(300,1e6) returns the derivative of the compressibility (at constant pressure) 
%                           at T =300 K and P = 1e6 Pa
%
%   6) H : Returns the enthalpy (J/kg) at a given Temperature (K) and Pressure (Pa)
%       ex: H(300,1e6) returns the enthalpy at T =300 K and P = 1e6 Pa
%   
%   7) Cv : Returns the specific heat at constant volume (J/kg/K) at a given Temperature (K) and Pressure (Pa)
%       ex: Cv(300,1e6) returns the enthalpy at T =300 K and Pa = 1e6 Pa
%
%   8)  Cp : Returns the specific heat at constant volume (J/kg/K) at a given Temperature (K) and Pressure (Pa)
%       ex: Cp(300,1e6) returns the enthalpy at T =300 K and Pa = 1e6 Pa
%
%   9) kH2 : Returns the conductivity (W/m/k) at a given Temperature (K) and Pressure (Pa)
%       ex: Cp(300,1e6) returns the enthalpy at T =300 K and Pa = 1e6 Pa
%
%   10) mu : Returns the viscosity (Pa*s) at a given Temperature (K) and Pressure (Pa)
%       ex: Cp(300,1e6) returns the enthalpy at T =300 K and Pa = 1e6 Pa
%
%   (we also return the R constant for hydrogen)
%
%
% ======================================================================
% Copyright (c) November 2023, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Build a scatter interpolant function for compressibility in the T;P
% domain

% Domain for interpolation: Tmin < T < Tmax 
%                           Pmin < P < Pmax

R = 4124.5; %J/kg/K

Tmin = -130; % (oC)
Tmax = 130;  % (oC)
Pmin = 1e4;  % (Pa)
Pmax = 100e6; % (Pa)
N = 100;      % (Number of data in each domain to build the Interpolation functions => N^2 Points)  

Z_TP = funcScatInterpol_Z_TP(Tmin,Tmax,Pmin,Pmax,N);

Z_Trho = funcScatInterpol_Z_Trho(Tmin,Tmax,Pmin,Pmax,N);

dZdT_rho = funcScatInterpol_dZdT_rho(Tmin,Tmax,Pmin,Pmax,N);

dZdT_P = funcScatInterpol_dZdT_P(Tmin,Tmax,Pmin,Pmax,N);

dZdP_T = funcScatInterpol_dZdP_T(Tmin,Tmax,Pmin,Pmax,N);

H = funcScatInterpol_H(Tmin,Tmax,Pmin,Pmax,N);

Cv = funcScatInterpol_Cv(Tmin,Tmax,Pmin,Pmax,N);

Cp = funcScatInterpol_Cp(Tmin,Tmax,Pmin,Pmax,N);

kH2 = funcScatInterpol_kH2(Tmin,Tmax,Pmin,Pmax,N);

mu = funcScatInterpol_mu(Tmin,Tmax,Pmin,Pmax,N);

save('ScatterInterpFunctions','R','Z_Trho','Z_TP','dZdT_rho','dZdP_T','dZdT_P','H','Cv','Cp','kH2','mu');

end

function Ho = funcEnthalpyLow(T)
    % https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1#Thermo-Gas
    % Chase, M.W., Jr., NIST-JANAF Themochemical Tables, Fourth Edition, J. Phys. Chem. Ref. Data, Monograph 9, 1998, 1-1951. [all data] 
    
    % T (K)
    % returns enthalpy (J/kg) at low pressure

    Href = 3931.7; %kJ/kg  % Enthalpy of hydrogen at 25 oC  (at low pressure)
    
    A = 33.066178;
    B = -11.363417;
    C = 11.432816;
    D = -2.772874;
    E = -0.158558;
    F = -9.980797;
    
    T = T/1000;
    
    Href = Href*2.01588e-3; % kJ/mol
    
    Ho = (Href + A*T + B*T.^2/2 + C*T.^3/3 + D*T.^4/4 - E/T + F)*1e3; %J/mol
    
    Ho = Ho/2.01588e-3; %(J/Kg)
end


function H = funcEnthalpy(T,P)
    % Eq 9 of Kushnir 2012 (Journal of Energy Resources)
    % The equation is modified in order to integrate in the Pressure
    % domain instead of the \rho domain 

    R = 4124.5;  %J/kg/K
    Ho = funcEnthalpyLow(T);
    %Z = funcZLemmon(T,P);
    
   % H = Ho + R.*T.*(Z-1)-integral(@(Pi) (R*T.^2.*funcdZdTLemmon(T,Pi)).*(1./Pi-1./funcZLemmon(T,Pi).*funcdZdPLemmon(T,Pi)),0,P);
    
    % see personnal note eq. 61 
    H = Ho - integral(@(Pi) (R*T.^2./Pi.*funcdZdT_P(T,Pi)),0,P);
end


%{
% I use instead the next function based on the link between Cp and Cv.
% It seems more stable to integrate cp instead of cv... 
% (don't know why ...?)

function Cv = funcCv(T,P)
    % Eq 10 of Kushnir 2012 (Journal of Energy Resources)
    % The equation is modified in order to integrate in the Pressure
    % domain instead of the \rho domain 
    % 
    % cv (J/Kg/K) 
        
     R = 4124.5; %(J/Kg/K);   
         
     A = 33.066178;
     B = -11.363417;
     C = 11.432816;
     D = -2.772874;
     E = -0.158558;
    
     % for Cp0 at T
     % https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1#Thermo-Gas
     % Chase, M.W., Jr., NIST-JANAF Themochemical Tables, Fourth Edition, J. Phys. Chem. Ref. Data, Monograph 9, 1998, 1-1951. [all data] 
     T = T/1000;
     Cp0 = A + B*T + C*T.^2 + D.*T^3 + E./T^2; %J/mol/K
     Cp0 = Cp0/2.01588e-3; %J/kg/K
     T = T*1000;

     Cv = Cp0 - R - R*T*integral(@(Pi) (T.*funcd2ZdTLemmon(T,Pi) + 2*funcdZdTLemmon(T,Pi)).*(1./Pi-1./funcZLemmon(T,Pi).*funcdZdPLemmon(T,Pi)),0,P);
       
end
%}


function cv = funcCv(T,P)

    Z = funcZLemmon(T,P);
    beta = 1./T + 1./Z.*funcdZdT_P(T,P);
    alpha = 1./P - 1./Z.*funcdZdP_T(T,P);
   
    R = 4124.5;
    cp = funcCp(T,P);
    
    cv = cp-(beta).^2.*T.^2.*Z*R./alpha./P;

end


function cp = funcCp(T,P)
   
    % Specific heat at constant pressure (see my personnal notes)
    % 
    % cp (J/Kg/K) 
        
     R = 4124.5; %(J/Kg/K);   
        
     cp0 = fcnCp0(T);

     cp = cp0 - R*T*integral(@(Pi) (T.*funcd2ZdT_P(T,Pi) + 2*funcdZdT_P(T,Pi)).*(1./Pi),0,P);
end


function Z = funcZLemmon(T,P)
    %   Returns compressibility factor
    %   for hydrogen according to Lemmon equation
    %
    %   P : Pressure (Pa)
    %   T : Temperature (K)
    %   Z : Compressibility factor
    %
    %
    
    N = 9;
    a = zeros(1,N);
    b = zeros(1,N);
    c = zeros(1,N);
    a(1) = 0.05888460;
    a(2) = -0.06136111;
    a(3) = -0.002650473;
    a(4) = 0.002731125;
    a(5) = 0.001802374;
    a(6) = -0.001150707;
 %   a(6) = -0.0012150707;
    a(7) = 0.958852e-4;
    a(8) = -0.1109040e-6;
    a(9) = 0.1264403e-9;
    
    b(1) = 1.325;
    b(2) = 1.87;
    b(3) = 2.5;
    b(4) = 2.8;
    b(5) = 2.938;
    b(6) = 3.14;
    b(7) = 3.37;
    b(8) = 3.75;
    b(9) = 4.0;
    
    c(1) = 1.0;
    c(2) = 1.0;
    c(3) = 2.0;
    c(4) = 2.0;
    c(5) = 2.42;
    c(6) = 2.63;
    c(7) = 3.0;
    c(8) = 4.0;
    c(9) = 5.0;
    Z = ones(size(P));

    P = P/1e6; 
    for i = 1:N
        Z = Z + a(i)*(100./T).^b(i).*P.^c(i);
    end

end



function dZdT = funcdZdT_rho(T,P)
        
    % Retruns dZ/dT at constant density as a function of T and P
    % 
    %   T: Temperature (K)
    %   P: Pressure (Pa)
    %   dZdT: The derivative of Z relative to T at constant rho

    N = 9;
    a = zeros(1,N);
    b = zeros(1,N);
    c = zeros(1,N);
    a(1) = 0.05888460;
    a(2) = -0.06136111;
    a(3) = -0.002650473;
    a(4) = 0.002731125;
    a(5) = 0.001802374;
    a(6) = -0.001150707;
 %   a(6) = -0.0012150707;
    a(7) = 0.958852e-4;
    a(8) = -0.1109040e-6;
    a(9) = 0.1264403e-9;

    b(1) = 1.325;
    b(2) = 1.87;
    b(3) = 2.5;
    b(4) = 2.8;
    b(5) = 2.938;
    b(6) = 3.14;
    b(7) = 3.37;
    b(8) = 3.75;
    b(9) = 4.0;
    
    c(1) = 1.0;
    c(2) = 1.0;
    c(3) = 2.0;
    c(4) = 2.0;
    c(5) = 2.42;
    c(6) = 2.63;
    c(7) = 3.0;
    c(8) = 4.0;
    c(9) = 5.0;

    
    Z = funcZLemmon(T,P);
    
    P = P/1e6;
    temp1 = zeros(size(P));
    temp2 = zeros(size(P));
    for i = 1:N
        temp = a(i)*(100)^b(i)./(T.^(b(i)+1)).*P.^c(i);
        temp1 = temp1 + b(i)*temp;
        temp2 = temp2 + c(i)*temp;
    end
    dZdT = (-temp1 + temp2)./(1-T.*temp2./Z);
end


function dZdT = funcdZdT_P(T,P)
        
    % Retruns dZ/dT at constant pressure as a function of T and P
    % 
    %   T: Temperature (K)
    %   P: Pressure (Pa)
    %   dZdT: The derivative of Z relative to T at constant Pressure

    N = 9;
    a = zeros(1,N);
    b = zeros(1,N);
    c = zeros(1,N);
    a(1) = 0.05888460;
    a(2) = -0.06136111;
    a(3) = -0.002650473;
    a(4) = 0.002731125;
    a(5) = 0.001802374;
    a(6) = -0.001150707;
 %   a(6) = -0.0012150707;
    a(7) = 0.958852e-4;
    a(8) = -0.1109040e-6;
    a(9) = 0.1264403e-9;

    b(1) = 1.325;
    b(2) = 1.87;
    b(3) = 2.5;
    b(4) = 2.8;
    b(5) = 2.938;
    b(6) = 3.14;
    b(7) = 3.37;
    b(8) = 3.75;
    b(9) = 4.0;
    
    c(1) = 1.0;
    c(2) = 1.0;
    c(3) = 2.0;
    c(4) = 2.0;
    c(5) = 2.42;
    c(6) = 2.63;
    c(7) = 3.0;
    c(8) = 4.0;
    c(9) = 5.0;

    dZdT = zeros(size(P));

    P = P/1e6;
    for i = 1:N
        dZdT = dZdT - b(i)*a(i)*100^b(i)./T.^(b(i)+1).*P.^c(i);
    end
end


function d2ZdT2 = funcd2ZdT_P(T,P)
        
    % Retruns d2Z/dT2 at constant pressure as a function of T and P
    % 
    %   T: Temperature (K)
    %   P: Pressure (Pa)
    %   d2ZdT2: The second derivative of Z relative to T at constant Pressure

    N = 9;
    a = zeros(1,N);
    b = zeros(1,N);
    c = zeros(1,N);
    a(1) = 0.05888460;
    a(2) = -0.06136111;
    a(3) = -0.002650473;
    a(4) = 0.002731125;
    a(5) = 0.001802374;
    a(6) = -0.001150707;
 %   a(6) = -0.0012150707;
    a(7) = 0.958852e-4;
    a(8) = -0.1109040e-6;
    a(9) = 0.1264403e-9;

    b(1) = 1.325;
    b(2) = 1.87;
    b(3) = 2.5;
    b(4) = 2.8;
    b(5) = 2.938;
    b(6) = 3.14;
    b(7) = 3.37;
    b(8) = 3.75;
    b(9) = 4.0;
    
    c(1) = 1.0;
    c(2) = 1.0;
    c(3) = 2.0;
    c(4) = 2.0;
    c(5) = 2.42;
    c(6) = 2.63;
    c(7) = 3.0;
    c(8) = 4.0;
    c(9) = 5.0;

    d2ZdT2 = zeros(size(P));

    P = P/1e6;
    for i = 1:N
        d2ZdT2 = d2ZdT2+ b(i)*(b(i)+1)*a(i)*100^b(i)./T.^(b(i)+2).*P.^c(i);
    end
end


function dZdP = funcdZdP_T(T,P)
    % Retruns dZ/dP at constant T from Lemmon equation as a function of T and P
    % 
    %   T: Temperature (K)
    %   P: Pressure (Pa)
    %   dZdP: The first derivative of Z relative to P 


    N = 9;
    a = zeros(1,N);
    b = zeros(1,N);
    c = zeros(1,N);
    a(1) = 0.05888460;
    a(2) = -0.06136111;
    a(3) = -0.002650473;
    a(4) = 0.002731125;
    a(5) = 0.001802374;
    a(6) = -0.001150707;
 %   a(6) = -0.0012150707;
    a(7) = 0.958852e-4;
    a(8) = -0.1109040e-6;
    a(9) = 0.1264403e-9;
    
    b(1) = 1.325;
    b(2) = 1.87;
    b(3) = 2.5;
    b(4) = 2.8;
    b(5) = 2.938;
    b(6) = 3.14;
    b(7) = 3.37;
    b(8) = 3.75;
    b(9) = 4.0;
    
    c(1) = 1.0;
    c(2) = 1.0;
    c(3) = 2.0;
    c(4) = 2.0;
    c(5) = 2.42;
    c(6) = 2.63;
    c(7) = 3.0;
    c(8) = 4.0;
    c(9) = 5.0;
    dZdP = zeros(size(P));

    P = P/1e6;
    for i = 1:N
        dZdP = dZdP + a(i)*c(i)*(100./T).^b(i).*P.^(c(i)-1);
    end
    dZdP = dZdP/1e6;
end

%{
function d2ZdT = funcd2ZdTLemmon(T,P)

% This function is not used anymore. (it was needed to calculate cv
% from integration. It seems to better more efficient to calculate
% cp from integration.
% 

    % Retruns d^2Z/dT^2 at constant density as a function of T and P
    % 
    %   T: Temperature (K)
    %   P: Pressure (Pa)
    %   d2ZdT: The second derivative of Z relative to T at constant rho


    N = 9;
    a = zeros(1,N);
    b = zeros(1,N);
    c = zeros(1,N);
    a(1) = 0.05888460;
    a(2) = -0.06136111;
    a(3) = -0.002650473;
    a(4) = 0.002731125;
    a(5) = 0.001802374;
    a(6) = -0.001150707;
 %   a(6) = -0.0012150707;
    a(7) = 0.958852e-4;
    a(8) = -0.1109040e-6;
    a(9) = 0.1264403e-9;

    b(1) = 1.325;
    b(2) = 1.87;
    b(3) = 2.5;
    b(4) = 2.8;
    b(5) = 2.938;
    b(6) = 3.14;
    b(7) = 3.37;
    b(8) = 3.75;
    b(9) = 4.0;
    
    c(1) = 1.0;
    c(2) = 1.0;
    c(3) = 2.0;
    c(4) = 2.0;
    c(5) = 2.42;
    c(6) = 2.63;
    c(7) = 3.0;
    c(8) = 4.0;
    c(9) = 5.0;

    
    dZdT = funcdZdTLemmon(T,P);
    Z = funcZLemmon(T,P);
    
    P = P/1e6;
    
    dPdT = P./T +P./Z.*dZdT;
    
    temp1 = zeros(size(P));
    temp2 = zeros(size(P));
    temp3 = zeros(size(P));
    temp4 = zeros(size(P));
    for i = 1:N
        temp1 = temp1 + a(i)*b(i)*(b(i)+1)*(100)^b(i)./(T.^(b(i)+2)).*P.^c(i);
        temp2 = temp2 + a(i)*b(i)*c(i)*(100)^b(i)./(T.^(b(i)+1)).*P.^(c(i)-1);
        temp3 = temp3 + a(i)*c(i)*(c(i)-1)*(100./T).^b(i).*P.^(c(i)-2);
        temp4 = temp4 + a(i)*c(i)*(100./T).^b(i).*P.^c(i);
    end
    d2ZdT = (temp1 - 2*temp2.*dPdT + temp3.*(dPdT).^2 + 2*temp4./Z./T.*dZdT)./(1-temp4./Z);
end
%}

function mu = funcViscosityMuzny(T,rho)
%Muzny et al 
% Correlation for the Viscosity of Normal Hydrogen Obtained from
% Symbolic Regression
% Journal of Chemical Engeneering data 2013, 58 pp 969-979
% dx.doi.org/10.1021/je301273j 

    a = zeros(1,5);
    b = zeros(1,7);
    c = zeros(1,6);
    a(1) = 2.09630e-1;
    a(2) = -4.55274e-1;
    a(3) = 1.43602e-1;
    a(4) = -3.35325e-2;
    a(5) = 2.76981e-3;
    b(1) = -0.1870;
    b(2) = 2.4871;
    b(3) = 3.7151;
    b(4) = -11.0972;
    b(5) = 9.0965;
    b(6) = -3.8292;
    b(7) = 0.5166;
    c(1) = 6.43449673;
    c(2) = 4.56334068e-2;
    c(3) = 2.32797868e-1;
    c(4) = 9.58326120e-1;
    c(5) = 1.27941189e-1;
    c(6) = 3.63576595e-1; 
    Na = 6.022137e23;
    
    M = 2.01588;
    sigma = 0.297;
    Ts = T/30.41;
    lnS = 0;
    for i = 1:5
        lnS = lnS + a(i)*log(Ts)^(i-1);
    end
    Ss = exp(lnS);
    
    eta0 = 0.021357*(M*T)^0.5/sigma^2/Ss; 
    
    Bs = 0;  
    for i = 1:7
        Bs = Bs + b(i)/(Ts)^(i-1);
    end
    B_eta = Na*Bs*(sigma*1e-9)^3;

    eta1 = B_eta*eta0;
    
    rho_sc = 90.909090909;
    T_c = 33.145;
    rho_r = rho/rho_sc;
    T_r= T/T_c;
    
    temp = c(2)*T_r+c(3)/T_r+c(4)*rho_r^2/(c(5)+T_r)+c(6)*rho_r^6; 
    mu = (eta0 + eta1*rho + c(1)*rho_r^2*exp(temp))*1e-6;

end



function lambda = funcThermalConductivityAssael(T,rho)
% Assael, M. , Assael, J. , Huber, M. , Perkins, R. and Takata, Y. (2011), 
% Correlation of the Thermal Conductivity of Normal and Parahydrogen from the 
% Triple Point to 1000 K and up to 100 MPa, 
% J. Phys. & Chem. Ref. Data (JPCRD), National Institute of Standards and Technology, 
% Gaithersburg, MD, [online], 
% https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=907928 (Accessed September 29, 2023) 
T_c = 33.145;
rho_c = 31.262;
A1 = zeros(1,7);
A2 = zeros(1,4);
B1 = zeros(1,5);
B2 = zeros(1,5);
A1(1) = -3.40976e-1;
A1(2) = 4.58820; 
A1(3) = -1.45080;
A1(4) = 3.26394e-1; 
A1(5) = 3.16939e-3;
A1(6) = 1.90592e-4;
A1(7) = 1.13900e-6;
A2(1) = 1.38497e2; 
A2(2) = -2.21878e1; 
A2(3) = 4.57151; 
A2(4) = 1;
B1(1) = 3.63081e-2;
B1(2) = -2.07629e-2; 
B1(3) = 3.14810e-2; 
B1(4) = -1.43097e-2; 
B1(5) = 1.74980e-3; 
B2(1) = 1.83370e-3;
B2(2) = -8.86716e-3;
B2(3) = 1.58260e-2;
B2(4) = -1.06283e-2;
B2(5) = 2.80673e-3;
C1 = 6.24e-4;
C2 = -2.58e-7;
C3 = 0.837;

temp1= 0;
temp2 =0;
T_r = T/T_c;
rho_r = rho/rho_c;
for i =1:7
    temp1 = temp1 + A1(i)*T_r^(i-1);
end
for i = 1:4
    temp2 = temp2 + A2(i)*T_r^(i-1);
end
temp = 0;
for i = 1:5
    temp = temp+ (B1(i)+B2(i)*T_r)*rho_r^i;
end
Delta_lambda_c = (C1/(C2+ abs(T_r-1)))*exp(-(C3*abs(rho_r-1))^2);


lambda0 = temp1/temp2;

lambda = (lambda0 + temp+Delta_lambda_c);
end



function Z_Trho = funcScatInterpol_Z_Trho(Tmin,Tmax,Pmin,Pmax,N)
    R = 4124.5;   
   
    Tmesh = linspace(Tmin,Tmax,N) + 273.15;
    Pmesh = linspace(Pmin,Pmax,N);
  %  Pmesh = logspace(2,8,500);
    [Tmesh,Pmesh] = meshgrid(Tmesh,Pmesh);
    Zmesh = funcZLemmon(Tmesh,Pmesh);
    
    rhomesh = Pmesh./Zmesh./Tmesh/R;   
    Zgrid = reshape(Zmesh,[],1);
    Tgrid = reshape(Tmesh,[],1);
    rhogrid = reshape(rhomesh,[],1);
    
    Z_Trho = scatteredInterpolant(Tgrid,rhogrid,Zgrid);
end

function dZdT = funcScatInterpol_dZdT_rho(Tmin,Tmax,Pmin,Pmax,N)

    Tmesh = linspace(Tmin,Tmax,N) + 273.15;
    Pmesh = linspace(Pmin,Pmax,N);
    [Tmesh,Pmesh] = meshgrid(Tmesh,Pmesh);
    dZdTmesh = funcdZdT_rho(Tmesh,Pmesh);
    dZdT = reshape(dZdTmesh,[],1);
    T = reshape(Tmesh,[],1);
    P = reshape(Pmesh,[],1);
    
    dZdT = scatteredInterpolant(T,P,dZdT);
end

function Z_TP = funcScatInterpol_Z_TP(Tmin,Tmax,Pmin,Pmax,N)
    %%%%%%%%%%%%%%%%%%
    %
    % Build a scatterInterpolant compressibility function using data from Lemmon equation
    %
    % Z(T,P)

    T = linspace(Tmin,Tmax,N) + 273.15;
    P = linspace(Pmin,Pmax,N);
    [T,P] = meshgrid(T,P);

    Z = funcZLemmon(T,P);
    
    Z = reshape(Z,[],1);
    T = reshape(T,[],1);
    P = reshape(P,[],1);

    Z_TP = scatteredInterpolant(T,P,Z);
end


function dZdT = funcScatInterpol_dZdT_P(Tmin,Tmax,Pmin,Pmax,N)

    Tmesh = linspace(Tmin,Tmax,N) + 273.15;
    Pmesh = linspace(Pmin,Pmax,N);
    [Tmesh,Pmesh] = meshgrid(Tmesh,Pmesh);
    dZdTmesh = funcdZdT_P(Tmesh,Pmesh); %dZdT at constant P
    dZdT = reshape(dZdTmesh,[],1);
    T = reshape(Tmesh,[],1);
    P = reshape(Pmesh,[],1);
    
    dZdT = scatteredInterpolant(T,P,dZdT);
end

function dZdP = funcScatInterpol_dZdP_T(Tmin,Tmax,Pmin,Pmax,N)

    Tmesh = linspace(Tmin,Tmax,N) + 273.15;
    Pmesh = linspace(Pmin,Pmax,N);
    [Tmesh,Pmesh] = meshgrid(Tmesh,Pmesh);
    dZdTmesh = funcdZdP_T(Tmesh,Pmesh); %dZdP at constant T
    dZdT = reshape(dZdTmesh,[],1);
    T = reshape(Tmesh,[],1);
    P = reshape(Pmesh,[],1);
    
    dZdP = scatteredInterpolant(T,P,dZdT);
end



function H = funcScatInterpol_H(Tmin,Tmax,Pmin,Pmax,N)
 
    T = linspace(Tmin,Tmax,N) + 273.15;
    P = linspace(Pmin,Pmax,N);
    [T,P] = meshgrid(T,P);
    T = reshape(T,[],1);
    P = reshape(P,[],1);
    
    N = length(T);
    H = zeros(N,1);
    for t = 1:N
        H(t) = funcEnthalpy(T(t),P(t));
    end
    
    H = scatteredInterpolant(T,P,H);
end

%{
function Cv = funcScatInterpol_Cv(Tmin,Tmax,Pmin,Pmax,N)
  
    T = linspace(Tmin,Tmax,N) + 273.15;
    P = linspace(Pmin,Pmax,N);
    [T,P] = meshgrid(T,P);
    T = reshape(T,[],1);
    P = reshape(P,[],1);
    N = length(T);
    
    Cv = zeros(N,1);
    for t = 1:N
        Cv(t) = funcCv(T(t),P(t));
    end
    
    Cv = scatteredInterpolant(T,P,Cv);
end
%}

function Cv = funcScatInterpol_Cv(Tmin,Tmax,Pmin,Pmax,N)
  
    T = linspace(Tmin,Tmax,N) + 273.15;
    P = linspace(Pmin,Pmax,N);
    [T,P] = meshgrid(T,P);
    T = reshape(T,[],1);
    P = reshape(P,[],1);
    N = length(T);
    
    Cv = zeros(N,1);
    for t = 1:N
        Cv(t) = funcCv(T(t),P(t));
    end
    
    Cv = scatteredInterpolant(T,P,Cv);
end



function Cp = funcScatInterpol_Cp(Tmin,Tmax,Pmin,Pmax,N)
  
    T = linspace(Tmin,Tmax,N) + 273.15;
    P = linspace(Pmin,Pmax,N);
    [T,P] = meshgrid(T,P);
    T = reshape(T,[],1);
    P = reshape(P,[],1);
    N = length(T);
    
    Cp = zeros(N,1);
    for t = 1:N
        Cp(t) = funcCp(T(t),P(t));
    end
    
    Cp = scatteredInterpolant(T,P,Cp);
end

function kH2 = funcScatInterpol_kH2(Tmin,Tmax,Pmin,Pmax,N)
% Interpolating function for heat conductivity of Hydrogen 
% based on Assel Regression model (see function : funcThermalConductivityAssael)

    R = 4124.5; 
    T = linspace(Tmin,Tmax,N) + 273.15;
    P = linspace(Pmin,Pmax,N);
    [T,P] = meshgrid(T,P);
    T = reshape(T,[],1);
    P = reshape(P,[],1);
    Z = funcZLemmon(T,P);
    rho = P./Z./T/R;
    N = length(T);
    kH2 = zeros(N,1);
    for t = 1:N
        kH2(t) = funcThermalConductivityAssael(T(t),rho(t));
    end
    kH2 = scatteredInterpolant(T,P,kH2);
end

function mu = funcScatInterpol_mu(Tmin,Tmax,Pmin,Pmax,N)
% Interpolating function for  viscosity of Hydrogen 
% based on Muzny Regression model (see function : funcViscosityMuzny(T,rho)
    R = 4124.5; 
    T = linspace(Tmin,Tmax,N) + 273.15;
    P = linspace(Pmin,Pmax,N);
    [T,P] = meshgrid(T,P);
    T = reshape(T,[],1);
    P = reshape(P,[],1);
    Z = funcZLemmon(T,P);
    rho = P./Z./T/R;
    N = length(T);
    mu = zeros(N,1);
    for t = 1:N
        mu(t) = funcViscosityMuzny(T(t),rho(t));
    end
    mu = scatteredInterpolant(T,P,mu);

end

