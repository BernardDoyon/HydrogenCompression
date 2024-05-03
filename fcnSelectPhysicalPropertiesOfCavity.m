function [flag_underground, flag_geometry, Rw_int, H,Dinlet,FunctThermalProperties] = fcnSelectPhysicalPropertiesOfCavity(flag)
    
   switch(flag)
       case 1 
           % Exampal 1 : Vertical tank with 3 different layers of different conductivities 
           flag_underground = 1;   % 0: Tank is outside (surrounded by air)  % 1: Tank is underground                                  
           flag_geometry = 2;      % 1: Horizontal tank     2: Vertical tank
           Rw_int = 0.25;          % radius (m)
           H = 200;                % height (m)
           Dinlet = 0.005;         % Inlet diameter for injection (m)
           FunctThermalProperties = @fcnThermalProperties3; %see bellow (function with 3 conductivities)
       case 2
           % Exampal 2 : Vertical tank with 2 different layers of different
           % conductivities (salt is one of the two laayers)
           flag_underground = 1;   % 0: Tank is outside (surrounded by air)  % 1: Tank is underground                                  
           flag_geometry = 2;      % 1: Horizontal tank     2: Vertical tank
           Rw_int = 0.25;          % radius (m)
           H = 200;                % height (m)
           Dinlet = 0.25;         % Inlet diameter for injection (m)
           FunctThermalProperties = @fcnThermalProperties2Salt;
       case 3
           % Exampal 3 : Vertical hole in salt
           flag_underground = 1;   % 0: Tank is outside (surrounded by air)  % 1: Tank is underground                                  
           flag_geometry = 2;      % 1: Horizontal tank     2: Vertical tank
           Rw_int = 0.25;          % radius (m)
           H = 200;                % height (m)
           Dinlet = 0.25;         % Inlet diameter for injection (m)
           FunctThermalProperties = @fcnThermalPropertiesSalt;
       case 4
           % Exampal 3 : Horizontal tank is outside (Couteau's article) 
           flag_underground = 0;   % 0: Tank is outside (surrounded by air)  % 1: Tank is underground                                  
           flag_geometry = 1;      % 1: Horizontal tank     2: Vertical tank
           Rw_int = 0.1294;          % radius (m)
           V = 0.036;
           H = V/pi/Rw_int^2;                % height (m)
           Dinlet = 0.003;         % Inlet diameter for injection (m)
           FunctThermalProperties = @fcnThermalPropertiesCouteau;          
   end
end






function [Rw, kRValue, CpValue, rhoValue]= fcnThermalProperties3(r,T)

  % T is not use!

    Rw = [0.35 .5 Inf]; % if only one thermal property => Rw = [Inf]
 
    kR(1) = 1.5;    % Thermal conductivity of first material from r*=1 (W/m/K)
    kR(2) = 2.5;    % Thermal conductivity of first material from r*=1 (W/m/K)
    kR(3) = 4;       % Thermal conductivity of first material from r*=1 (W/m/K)
    
    CpR(1) = 1.5e3;       % Specific heat of first material (J/Kg/K)
    CpR(2) = 0.84e3;      % Specific heat J/Kg/K of second material (J/Kg/K)
    CpR(3) = 1.2e3;      % Specific heat J/Kg/K of second material (J/Kg/K)
   
    rhoR(1) = 2.6e3;      % Density of first material (Kg/m^3)
    rhoR(2) = 2e3;        % Density of second material (Kg/m^3)
    rhoR(3) = 3e3;        % Density of second material (Kg/m^3)
    
    if r<=Rw(1)
        kRValue = kR(1);
        CpValue = CpR(1);
        rhoValue = rhoR(1);
    elseif r<=Rw(2)
        kRValue =  kR(2);
        CpValue = CpR(2);
        rhoValue = rhoR(2);
    else
        kRValue =  kR(3);
        CpValue = CpR(3);
        rhoValue = rhoR(3);
    end

end




function [Rw, kRValue, CpValue, rhoValue]= fcnThermalProperties2Salt(r,T)

    Rw = [0.35 Inf]; % if only one thermal property => Rw = [Inf]
 
    kR(1) = 1.5;      % Thermal conductivity of first material from r*=1 (W/m/K)
  %  kR(2) = variable with T; % don't need to set because second medium is salt and 
                              % we use a function of T instead (see bellow)
    
    CpR(1) = 1.5e3;       % Specific heat of first material (J/Kg/K)
    CpR(2) = 0.84e3;      % Specific heat J/Kg/K of second material (J/Kg/K)
   
    rhoR(1) = 2.6e3;      % Density of first material (Kg/m^3)
    rhoR(2) = 2e3;        % Density of second material (Kg/m^3)
    
    if r<=Rw(1)
        kRValue = kR(1);
        CpValue = CpR(1);
        rhoValue = rhoR(1);
    elseif r<=Rw(2)
        kRValue = fcnkRSalt(T);
        CpValue = CpR(2);
        rhoValue = rhoR(2);
    end

end

function [Rw, kRValue, CpValue, rhoValue]= fcnThermalPropertiesCouteau(r,T)

    % T is not use!

    Rw = [0.1354 0.1604];  % set the R values where the thermal properties are changing 
                           % if only one thermal property => Rw = [Inf]
 
    kR(1) = 0.36;    % Thermal conductivity of first material (from r*=1) (W/m/K)
    kR(2) = 1.5;     % Thermal conductivity of second material (from r*=1) (W/m/K)
    CpR(1) = 1880;   % Specific heat of first material (J/Kg/K)
    CpR(2) = 1400;   % Specific heat of second material (J/Kg/K)
    rhoR(1) = 947;      % Density of first material (Kg/m^3)
    rhoR(2) = 1600;
    
    if r<=Rw(1)
        kRValue = kR(1);
        CpValue = CpR(1);
        rhoValue = rhoR(1);
    else
        kRValue =  kR(2);
        CpValue = CpR(2);
        rhoValue = rhoR(2);
    end
end


function [Rw, kRValue, CpValue, rhoValue]= fcnThermalPropertiesSalt(r,T)

    Rw = [Inf]; % if only one thermal property => Rw = [Inf]

  %  kR(1) = variable with T;     % Thermal conductivity of first material from r*=1 (W/m/K)
    
    CpR(1) = 0.84e3;      % Specific heat J/Kg/K of second material (J/Kg/K)
   
    rhoR(1) = 2e3;        % Density of second material (Kg/m^3)
     
    
    if r<=Rw(1)
        kRValue = fcnkRSalt(T);
        CpValue = CpR(1);
        rhoValue = rhoR(1);
    end
    
end

function kR = fcnkRSalt(T)
    kR = -0.835*log(T)+8.732;
end