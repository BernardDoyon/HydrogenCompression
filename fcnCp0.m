function cp0 = fcnCp0(T)

    % for Cp0 at T
     % https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1#Thermo-Gas
     % Chase, M.W., Jr., NIST-JANAF Themochemical Tables, Fourth Edition, J. Phys. Chem. Ref. Data, Monograph 9, 1998, 1-1951. [all data] 
         
     A = 33.066178;
     B = -11.363417;
     C = 11.432816;
     D = -2.772874;
     E = -0.158558;
    
     % for Cp0 : 
     T = T/1000;
     cp0 = A + B*T + C*T.^2 + D.*T^3 + E./T^2; %J/mol/K
     cp0 = cp0/2.01588e-3; %J/kg/K
    
end