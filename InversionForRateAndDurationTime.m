%
%   This file calculates the injection rate and duration time needed for a
%   a target total mass to be injected in order to respect a constraint on
%   the maximal temperature reached by the gaz during filling. 
%
%   The constraint on the max temperature is set at line 53
%   The total mass to be injected is set at line 54
%
%   We first define a 2D function gx
%
%           gx(1) = Tmax(x1,x2) - TargetMaxTemperature
%           gx(2) = TotalMass(x1,x2) - TargetTotalMass
%
%           where   x1:     Injection duration time (s)
%                   x2:     Injection rate (kg/s)
%                   Tmax(x1,x2): The simulated maximal temperature with
%                                    x1 and x2
%                   TotalMass(x1,x2): The total mass injected (x1*x2)                                      
%
%   The solution corresponds to the zeros of gx. 
%   We use a damped Newton-Raphson iteration algorithm to find the zeros
%       (see ref: Davidchack and Lai,  Phys. Rev. E 60, 6172 – November 1999)
%       It was proven that this particular damped
%       correction is converging (slowly) to the solution when initial guess 
%       is far from it. When close to the solution, the damped correction
%       is almost equivalent to a Newton-Raphson step (here a 2D Newton-Raphson step)
%   In order to accelerate the convergence, we impose to use a 2D
%   Newton-Raphson step when the solution is within 5% of the Target
%   values. 
%     
%   The function gx can be calculated with:
%
%       1- A simulated injection with real hydrogen gaz properties
%       2- A simulated injection with simple hydrogen gaz properties 
%               (Z = 1 or Z = Z_T0)
%       3- An estimated max temperature from the analytical solution (Perfect conductor is assumed)  
%
%         (we use a function handle @ for gx) 
%
%
%   The execution of the file produces a graph that shows the solutions
%   (for real gaz, for simple gaz and with analytical expression when perfect
%   conductir is assumed). The solutions are marked with a symbol  
%   on the target mass contour plot.
%
%
% ======================================================================
% Copyright (c) October 2022, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================

Traget = zeros(2,1);
x = zeros(2,1);
Target(1) = 70;   % The constraint: Maximal temperature during injection 
Target(2) = 1200; % Total mass to be injected
x(1,1) = 1*3600; %6*3600; %inital guess for injection duration time
x(2,1) = 0.02; %initial guess for injection rate


% *********** SET NUMERICAL PARAMETERS FOR SIMULATION ********************
%
%   Set all the numerical parameterrs and numerical solver options for the simulation:
%
%       - The cavity dimension parameters (radius and height)
%       - The gas parameters
%       - The rock parameters
%       - The Numerical solver options
%
for n=1:1
    % The cavity dimensions parameters

    p.underground = 1;            % 0: Tank is outside (surrounded by air)    
                                  %  1: Tank is underground
                                  
    p.geometry = 2;               % 1: Horizontal tank     2: Vertical tank
    
    p.Rw_int = 0.25;                 % radius (m)
    p.H = 200;                    % height (m)
    
    p.Ac_int = 2*pi*p.Rw_int*p.H;        % Surface of cylindrical cavity (m^2)
    p.V = pi*p.Rw_int^2*p.H;  
    
    p.Dinlet = 0.003;               % Inlet diameter for injection (m)
  
    p.model = 1;                % 0 -> Ideal gaz: Z = 1
    % 1 -> Simplified model: Z = Z0; dZdT = ZT0

    % The gas parameters

    p.R = 4124.5;           % Specific gas constant for hydrogen (J/Kg/K)
    p.Ti = 70 + 273.15;     % Injected gas temperature (K)
    p.T0 = 30 + 273.15;     % Initial gas temperature (K)
    p.P0 = 1.0133e5;        % Initial gas pressure (Pa)
    p.v = 0;                % Gas speed in inlet when charging and discharging (m/s)
    p.dz = p.H/2;           % Vertical distance between inlet and middle of the tank (m)
    
    % The rock parameters

     p.TRw0 = p.T0;        % Initial rock temparature (K)
    %p.Pe = p.P0;          % Air pressur at reservoir edge (Pa)
    %p.rhoR = 2100;        % Rock density (kg/m^3)
    %p.cpR = 840;          % Constant pressure specific heat of rock (J/kg/K)
    p.hc_int = 30;         % Heat transfer coefficient  (W/m2/K)
    
    % Conductivity of material in tank from p.Rw_in to Inf
    p.kR(1) = 1;  %        % thermal conductivity of first material from r*=1 (W/m/K)
    %p.kR(2) = 2.4;         % Thermal conductivity of second material (W/m/K)
        
    
    p.alpha(1) = p.kR/2100/840;  % Thermal diffusivity of first material (m^2/s)
    % p.alpha(2) = 9.27e-7;   % Thermal diffusivity of second material (m^2/s)
    
    p.Rw = [Inf]; % if only one thermal property => p.Rw = [Inf]
    % The radius value where the material is changing in the tank
    %  1) from p.Rw_in to p.Rw(1) -> p.rho(1); p.cp(1); p.kR(1)
    %  2) from p.Rw(1) to p.Rw(2) -> p.rho(2); p.cp(2); p.kR(2)
    %               :::
    %  3) from p.Rw(N-1) to Inf -> p.rho(N); p.cp(N); p.kR(N)
    
    % In this comment, N is the number of different material in the tank

    
    % Numerical solver options

    SolverOptions.solver = 1;       % 0- ode45 solver (Non-stiff explicit Runge Kutta)
                                    % 1- ode23s solver (stiff solver)

    SolverOptions.N = 20;           % Number of grid point for the numerical solution in the rock
    SolverOptions.beta = 1.05;      % Beta parameter for the variable change r -> eta
    SolverOptions.Rp = 0;           % The minimum value of R in the rock where the temperature is not affected by heat conduction
                                    % When set to 0, eq A.7 of Kushnir 2012 is used
    SolverOptions.eps = 0.005;      % Parameter epsilon in eq A.7 of Kushnir
    SolverOptions.RelTol = 1e-4;    % Relative Tolerance for the numerical solver
    SolverOptions.AbsTol = 1e-6;    % Absolute Tolerance for the numerical solver
    SolverOptions.NtimeStep = 0;    % If NtimeStep = 0, the solution is return at every time step of the ode solver 
                                    % If NtimeStep = N, the solution is return for N time step of equal values.  
end
%
%% ********* END SET NUMERICAL PARAMETERS FOR SIMULATION ******************

%% ********************** SET INJECTION SCENARIO HERE *********************
%
for n= 1:1
 % The injection scenario parameters

    tp = linspace(0.15,6,10)*3600; %duration time scenarios (seconds)
    p.mc = 0;        % Mass rate variation during injection/extraction (kg/s) (for dm/dt simulations)
    p.pc = 0;     % Pressure rate variation during injection/extraction (Pa/s) (for dP/dt simulations)


    %   See file InjectStoraExtractScenarioFunction.m
    %
    %   Here, we need to pass a function handler to the variables p.Fie and
    %   p.Gie:
    %           p.Fie : Is the function handler that will return the value of
    %                   Fie (the modulation of the mass injection/extraction rate) for the time t*.
    %
    %           p.Gie : Is the function handler that will return the value of
    %                   Gie (the modulation of the pressure variation rate
    %                   during extraction/injection) for the time t*
    %
    %   Different functions are already defined in the file InjectStoraExtractScenarioFunction.m
    %   and this file can be modified in order to adapt any
    %   injection/stroage/extraction scenarios. We select the function handlers
    %   for p.Fie and p.Gie by passing two values to
    %   "InjectStoraExtractScenarioFunction"
    %
    %   Fie = 1 => Constant mass injection rate
    %   Fie = 2 => Constant mass injection rate from t* = 0 to t*= 1/2 and storage from t* = 1/2 to t* = 1
    %   Fie = 3 => Injection/Storage/Extraction/Strorage scneario
    %   Fie = 4 => Injection only but with mass injection rate decreasing linearly
    %
    %   Gie = 1 => Injection only from t* = 0 to t* = 1 with pressure ramp (dP/dt = constant)
    %   Gie = 2 => Injection from t*=0 to t*=1/2 (with dP/dt = constant) and storage from t* = 1/2 to t*= 1
    %   Gie = 3 => Injection/Storage/Extraction/Strorage scneario
    %
    %

    NbCycles = 1;
    Fie = 1;
    Gie = 1; 
end
%
%% ************************ END OF INJECTIOM SCENARIO *********************



%% *** Main 
%

if ~isfile('ScatterInterpFunctions.mat')
    fprintf('Generating interpolation functions and saving to ScatterInterpFunctions.mat \n');
    fcnGenerateHydrogenPropertiesFunctions()
end

[p.Fie, p.Gie] = fcnInjectStoraExtractScenarioFunctions(Fie,Gie); % function handlers for p.Fie and p.Gie

solution = zeros(3,4); % to store the solution
        % i = 1 Complex real gaz model
        % i = 2 Simple or ideal gaz model (depending on p.model)
        % i = 3 Analytical solution (perfect conductivity kr = \infty is assumed)
    
     % solution(i,1) is the injection duration rate (s) found by Newton-Raphson
     % solution(i,2) is the injection rate (kg/s) found by Newton-Raphson to respect the constraint on the max temperature
     % solution(i,3) is the max temperature with injection time duration and rate

 % The constraints (objective function)
 %  1) Maximal temperature during injection
 %  2) The total mass 

eps = 1e-1; % convergence criteria for the norm of gx 
maxite = 20;% maximum number of iteration

for i = 1:3 
    switch i
        case 3
            flagRealGas = 1;
            gxfunction = @(x) gxfunctionNumericalGasModel(x,Target,p,SolverOptions,flagRealGas);
            fprintf('Solution for the simulated temperature of Real gaz\r');
        case 2
            flagRealGas = 0;
            gxfunction = @(x) gxfunctionNumericalGasModel(x,Target,p,SolverOptions,flagRealGas);
            fprintf('Solution for the simulated temperature of simple gaz \r');
        case 1
            gxfunction = @(x) gxfunctionAnalytical(x,Target,p);
            fprintf('Solution for the max temperature estimated from analytical expression (Perfect conductor) \r');
    end
    
    gx(1) = 1;      % gx is the value return for the gx function for
    gx(2) = 1;      % which a zero is searched
                    % We set it to 1 at the beginning of the inversion process
   
    absgx = 1;      % absgx is the norm of gx: absgx = sqrt(g(x)'*g(x))
                    % We set it to 1 at the beginning of the inversion process
    iterate = 0;    % the counter for the number of iteration
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %    In the next while loop, we try to find the injection rate for the
    % proposed injection duration time in order to respect the constraint
    % on the maximal temperature
    %
    %  We use a danmped Newton-Raphson correction dx 
    %
    %   dx =  (ones(2,2)*beta*sqrt(obj)-dgdx)\gx;
    %
    %   with  beta ~ 0.0005
    %        
    %     - When close to the solution, this correction is equivalent to
    %        a 2D Newton-Raphson.
    %     - When far from the solution, it was proven that this particular damped
    %     correction is converging (slowly) to the solution  
    %       see ref: Davidchack,  Phys. Rev. E 60, 6172 – Published 1
    %       November 1999)
    %
    %    To accelerate the convergence, we force a Newton-Raphson correction if 
    %       the solution gx is within 5% of the Target values
    %
    %
    
    flag = 0;
    fprintf('Iteration for inversion process = 1');
    beta = 0.0005;
    while((absgx > eps)&&(iterate < maxite))
        iterate = iterate +1;
        if (iterate>1)
            for j=0:log10(iterate-1)
                fprintf('\b'); % delete previous counter display
            end
            fprintf('%d', iterate);
        end
        gx = gxfunction(x); %2D gx function
        absgx = sqrt((gx)'*gx); % norm of gx
        dgdx = dgdxfunction(gxfunction,x,gx); %jacobian matrix of gx
        if (abs(gx(1))>0.05*Target(1))||(abs(gx(2))>0.05*Target(2))
            dx = (ones(2,2)*beta*absgx-dgdx)\gx; %2D Damped Newton-Raphson 
        else
            dx = -dgdx\gx; %We are close to solution => 2D Newton-Raphson 
        end
        x = x + dx;
    end
    solution(i,1) = x(1);      %solution for injection duration time 
    solution(i,2) = x(2);      %solution for injection rate 
    solution(i,3) = gx(1) + Target(1); %The max temperature during injection for the calculated injection rate
    solution(i,4) = x(1)*x(2); %The injected mass found
    solution(i,5) = iterate; % number of iterates for the Newton-Raphson procedure
    fprintf('\r');
end
fprintf('\r');

%%%%%%%%%%%%%%%%
%
% Plotting the results with somme contour plots for total mass
%
%
dtAnalytical = solution(1,1);
dtSimplegaz = solution(2,1);
dtRealgaz = solution(3,1);

mcAnalytical = solution(1,2);
mcSimplegaz = solution(2,2);
mcRealgaz = solution(3,2);

TAnalytical = solution(1,3); 
TSimplegaz = solution(2,3); 
TRealgaz = solution(3,3); 


plot(dtRealgaz/3600, mcRealgaz,'*'); 
xlabel('Injection time (h)');
ylabel('${\dot m_c}$ (kg/s)','interpreter', 'latex','FontSize',16);
%ylim([0 0.3]);
hold on

if (p.model)
    textLegend = '(Simple gaz model)';
else
    textLegend = '(Ideal gaz model)';
end
plot(dtSimplegaz/3600, mcSimplegaz,'x'); 
plot(dtAnalytical/3600, mcAnalytical,'o'); 


%choose y limit according to analytical scenario
ylim([0 round(max(mcAnalytical),1)+0.1]);

dt = linspace(.15, 12, 100);
mc = Target(2)./(dt*3600);

dt = linspace(.15, 12, 100);

plot(dt,mc,'--');

legend(strcat('T_{Max} = ',num2str(TRealgaz),'^oC (Real gaz)'),...
    strcat('T_{Max} = ',num2str(TSimplegaz),strcat('^oC',textLegend)),...
    strcat('T_{Max} = ',num2str(TAnalytical),'^oC (Perfect conduction)'),...
    strcat('M =',num2str(Target(2)),'kg'));
hold off



function gx = gxfunctionNumericalGasModel(x,Target,p,SolverOptions,flagRealGas) %#ok<DEFNU>
    % Calculate gx with a Real gaz model
    %
    % Input
    %  1)Target (a vector):   Target(1) = Target temperature
    %                        Target(2) = Target Mass
    %  2) x (a vector): x(1) = Injection time duration (s)
    %                   x(2) = injection rate (kg/s)
    % 
    % Output:
    %   gx is a vector with 2 components
    %
    % gx(1) = Tmax - Target(1)
    % gx(2) = M - Target(2)
    % 
    % where 
    %   - Target(1) is the constraint on the maximal temperature
    %   - Tmax is the maximal temperature simulated during a injection
    %   scenario, the simulation is done with a rela gaz model 
    %   (using fcnSolveKushnirNum.m)
    %   - M is the total mass injected 
    %   - Target(2) is the target mass


    NbCycles = 1;
    p.tp = x(1);
    p.mc = x(2);
    para = fcnParametersForSolver(p);
    [~,~, Tair, ~, ~,~,~] = fcnSolve_dmdt(NbCycles,para,SolverOptions,flagRealGas);
    gx(1,1) = max(Tair)*para.T0-273.15 - Target(1);
    gx(2,1) = x(1)*x(2) - Target(2);
end




function gx = gxfunctionAnalytical(x,Target,p)

% Calculate gx with analytical solution

 % Input
    %  1)Target (a vector):   Target(1) = Target temperature
    %                        Target(2) = Target Mass
    %  2) x (a vector): x(1) = Injection time duration (s)
    %                   x(2) = injection rate (kg/s)
    % 
    % Output:
    %   gx is a vector with 2 components
    %
    % gx(1) = Tmax - Target(1)
    % gx(2) = M - Target(2)
    % 
    % where 
    %   - Target is the constraint on the maximal temperature
    %   - Tmax is the maximal temperature of hydrogen at the end of injection 
    %       based on analytical solution 
    %       (perfect conductor is assume 
    %       equation B.1 of Kushnir 2012 (Vol 55 Journal of Heat and mass transfer) 
    %    - M is the total mass injected 
    %    - Target(2) is the target mass
    
    p.td = x(1);
    p.tp = x(1);
    p.mc = x(2);

    para = fcnParametersForSolver(p);

    T0 = para.T0;
    TR = para.TRw0*T0;
    mr = para.mr;
    Ti = para.Ti*T0;
    gamma = para.gamma;
    qr = para.qr/para.cv0;  
    Rs = para.Rs;
    Us = para.Us;
    
    a1 = (gamma*Ti/TR + qr)/(gamma-Rs+qr);
    a2 = Us/(gamma-Rs+qr+1);
    b1= gamma-Rs+qr;
    
    rho = 1+mr;
    
    gx(1,1) = (a1+a2.*rho + (T0/TR - a1-a2)*rho.^(-b1))*TR - 273.15 - Target(1);
    gx(2,1) =  x(1)*x(2) - Target(2);
end 


function dgdx = dgdxfunction(gxfunction,x,gx)
    % This function calculates de Jacobian matrix of gx
    %
    % Jacobian matrix =  dg1/dx1  dg1/dx2
    %                    dg2/dx1 dg2/dx2
    %
    % (first line is evaluated numericcally)
    % (second line is obtained analytically)
    
    % x(1): Injection duration time
    % x(2): Injection rate

    dgdx = zeros(2,2);
    eps = 0.1;
    delta = 1e-3;
    for i = 1:2
        xp = x;
        dx = x(i)*eps+delta;
        xp(i) = x(i) + dx;
        gxp = gxfunction(xp);
        dgdx(1,i) = (gxp(1)-gx(1))/dx;
    end
    dgdx(2,1) = x(2);
    dgdx(2,2) = x(1);
end



