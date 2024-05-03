%
%   This file calculates the injection rate needed for a particular
%   injection duration time scenraio in order to respect a constraint on
%   the maximal temperature reach by the gaz during filling.
%
%  The constraint on the max temperature is set at line 42
%
%   We first define a function "gxfunction(x)"
%
%           gxfunction(x) = Tmax(x) - TargetMaxTemperature
%        
%           where   x:     Injection rate (kg/s)
%                   Tmax(x): The simulated maximal temperature for
%                             injection rate at x                                     
%
%   The solution corresponds to the zeros of gxfunction => gxfunction(x_solution) = 0 
%   We use a Newton-Raphson iteration algorithm to find the zero.
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
%   The execution of the file produces a graph showing the injection rate 
%   for 10 injection duration times (a star on the graph) respecting the
%   contraint in max temperature
%
%   Superposed on the graph are the contour plot graph for the total mass
%   injected regardless of the temperature.
%
%
% ======================================================================
% Copyright (c) October 2022, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================


Ttarget = 70;   % The constraint: Maximal temperature during injection 

%% *********** SET NUMERICAL PARAMETERS FOR SIMULATION ********************
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

    p.Rw_int = 0.25;              % radius (m)
    p.H = 200;                    % height (m)
    
    p.Ac_int = 2*pi*p.Rw_int*p.H;        % Surface of cylindrical cavity (m^2)
    p.V = pi*p.Rw_int^2*p.H;  
    
    p.Dinlet = 0.003;               % Inlet diameter for injection (m)
  
    p.model = 1;                % 0 -> Ideal gaz: Z = 1
    % 1 -> Simplified model: Z = Z0; dZdT = ZT0

    % The gas parameters

    p.R = 4124.5;           % Specific gas constant for hydrogen (J/Kg/K)
    p.Ti = 30 + 273.15;     % Injected gas temperature (K)
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

    
    % p.kr; p.phi and p.mu are not used yet in the code...
   % p.kr = 0;%5e-14;      % Calibrated permeability of rock (m^2) SET TO ZERO FOR NO LEAK SCENARIO
   % p.phi = 0.1;          % Porosity
   % p.mu = 1.79e-5;       % Air viscosity (Pa s)

    
    % Numerical solver options

    SolverOptions.solver = 1;       % 0- ode45 solver (Non-stiff explicit Runge Kutta)
                                    % 1- ode23s solver (stiff solver)

    SolverOptions.N = 140;           % Number of grid point for the numerical solution in the rock
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


%% ******** TASKS TO DO FOR THIS FILE ********************************
% 
% SOLVE_EQ -> do the numerical iteration 
% FIGURES - Plot options (to show the results) 
%               The figures are saves in "Figures" directory. 
%               This directory will be created in the parent directory
%
SOLVE_EQ = 1; %  %(0-> no; 1->yes)

% Figure 1 shows: the solution for the injection rate has a function of
%           the injection time in order to respect the limit temperature

FIGURE1 = 1; %(0->no; 1->yes)
Figure1Title = 'Injection rate for different injection time';
Figure1FileName = 'Figure1';
%
%%% End TASKS

%%
if(SOLVE_EQ) % The numerical solutions are evaluated here
    for n= 1:1 %#ok<UNRCH>
        if ~isfile('ScatterInterpFunctions.mat')
            fprintf('Generating interpolation functions and saving to ScatterInterpFunctions.mat \n');
            fcnGenerateHydrogenPropertiesFunctions()
        end
        
        [p.Fie, p.Gie] = fcnInjectStoraExtractScenarioFunctions(Fie,Gie); % function handlers for p.Fie and p.Gie
        
        solution = zeros(3,3,length(tp)); % to store the solution
        % i = 1 Complex real gaz model
        % i = 2 Simple or ideal gaz model (depending on p.model)
        % i = 3 Analytical solution (perfect conductivity kr = \infty is assumed)
        
        % solution(i,1,:) is the injection rate found by Newton-Raphson to respect the constraint on the max temperature
        % solution(i,2,:) is the max temperature with injection time duration and rate
        % solution(i,3,:) is the number of iterates for the Newton-Raphson procedure to converge
        
        epsAbs = 1e-1; % absolute convergence criteria for the objective function (temperature)
        epsRel = 1e-3; % rlelative precision on injection rate
        
        maxite = 10; % maximum number of iteration
        
        x = 0.02; %initial guess for the first injection scenario td = 0.15*3600 s
        
        for i = 1:3
            
            for j = 1:length(tp)
                fprintf('\r Injection scenario = %d; ',j);
                p.tp = tp(j);
                switch i
                    case 1
                        gxfunction = @(x) gxfunctionNumericalGasModel(x,Ttarget,p,SolverOptions,1);
                        fprintf('Solution for Real gaz scenarios \r');
                    case 2
                        gxfunction = @(x) gxfunctionNumericalGasModel(x,Ttarget,p,SolverOptions,0);
                        fprintf('Solution for simple gaz scenarios \r');
                    case 3
                        gxfunction = @(x) gxfunctionAnalytical(x,Ttarget,p);
                        fprintf('Solution for analytical temperature scenarios \r');
                end
                
                gx = 1;         % gx is the value return be the gx function
                % We set it to 1 at the beginning of the inversion process
                iterate = 0;    % the counter for the number of iteration
                
                % we try to find with the next while loop the injection rate for the
                % proposed injection duration time in order to respect the constraint
                % on the maximal temperature
                
                % Simple Newton-Raphson procedure
                fprintf('Iteration for inversion process = 1');
                dx = 1;
                while((abs(gx) > epsAbs)&&(iterate < maxite)&&(abs(dx/x)>epsRel))
                    iterate = iterate +1;
                    if (iterate>1)
                        for n=0:log10(iterate-1)
                            fprintf('\b'); % delete previous counter display
                        end
                        fprintf('%d', iterate);
                    end
                    gx = gxfunction(x);
                    dgdx = dgdxfunction(gxfunction,x,gx);
                    dx = -gx/dgdx;  % Newton-Raphson correction step (one-dimension)
                    x = x + dx;     % next iterate
                end
                solution(i,2,j) = x;      %solution for injection rate (will be the initial guess for next injection scenario)
                solution(i,3,j) = gx + Ttarget; %The max temperature during injection for the calculated injection rate
                solution(i,4,j) = iterate; % number of iterates for the Newton-Raphson procedure
            end
        end
    end
end
fprintf('\r');
%
%%

%%
% Plotting the results with somme contour plots for total mass
%
%
if (FIGURE1)
    for n= 1:1
        currentFolder = pwd;
        FigurePath = 'Figures';
        if ~isfolder(FigurePath)
            mkdir(FigurePath)
        end
        mcRealgaz = reshape(solution(1,2,:),[1,length(tp)]);
        mcSimplegaz = reshape(solution(2,2,:),[1,length(tp)]);
        mcAnalytical = reshape(solution(3,2,:),[1,length(tp)]);
        fig1= figure(1);
        plot(tp/3600, mcRealgaz,'-*');
        xlabel('Injection time (h)');
        ylabel('${\dot m_c}$ (kg/s)','interpreter', 'latex','FontSize',16);
        %ylim([0 0.3]);
        hold on
        
        if (p.model)
            textLegend = '(Simple gaz model)';
        else
            textLegend = '(Ideal gaz model)';
        end
        plot(tp/3600, mcSimplegaz,'-x');
        plot(tp/3600, mcAnalytical,'-o');
        
        
        %choose y limit according to analytical scenario
        ylim([0 round(max(mcAnalytical),1)+0.1]);
        
        
        M = linspace(200, 1400,7); % contour plots for total Mass
        dt = linspace(.15, 6, 100);
        for m = 1:length(M)
            mc = M(m)./(dt*3600);
            if (m == 1)
                plot(dt,mc,'--');
            else
                plot(dt,mc,'--');
            end
        end
        legend(strcat('T_{Max} = ',num2str(Ttarget),'^oC (Real gaz)'),...
            strcat('T_{Max} = ',num2str(Ttarget),strcat('^oC',textLegend)),...
            strcat('T_{Max} = ',num2str(Ttarget),'^oC (Perfect conduction)'),...
            strcat('M =',num2str(M(1)),'kg'),...
            strcat('M =',num2str(M(2)),'kg'),...
            strcat('M =',num2str(M(3)),'kg'),...
            strcat('M =',num2str(M(4)),'kg'),...
            strcat('M =',num2str(M(5)),'kg'),...
            strcat('M =',num2str(M(6)),'kg'),...
            strcat('M =',num2str(M(7)),'kg'));
        hold off
    end
    delete (strcat(currentFolder,'\',FigurePath,'\',Figure1FileName,'.png'));
    print(fig1,'-dpng', strcat(currentFolder,'\',FigurePath,'\',Figure1FileName));
end
%
%%

function gx = gxfunctionNumericalGasModel(x,Ttarget,p,SolverOptions,flagRealGas) %#ok<DEFNU>
    % Calculate gx with a Real gaz model
    %
    % gx = Tmax - Ttarget
    % 
    % where 
    %   - Ttarget is the constraint on the maximal temperature
    %   - Tmax is the maximal temperature simulated during a injection
    %   scenario, the simulation is done with a real gaz model 
    %   (using fcnSolveKushnirNum.m)

    NbCycles = 1;
    p.mc = x;
    para = fcnParametersForSolver(p);
    [~,~, Tair, ~, ~,~,~] = fcnSolve_dmdt(NbCycles,para,SolverOptions,flagRealGas);
    gx = max(Tair)*para.T0-273.15 - Ttarget;
end


function gx = gxfunctionAnalytical(x,Ttarget,p) %#ok<DEFNU>

% Calculate gx with analytical solution
    %
    % gx = Tmax - Ttarget
    % 
    % where 
    %   - Ttarget is the constraint on the maximal temperature
    %   - Tmax is the maximal temperature of hydrogen at the end of injection 
    %       based on analytical solution 
    %       (perfect conductor is assumed 
    %       equation B.1 of Kushnir 2012 (Vol 55 Journal of Heat and mass transfer) 
    
    p.mc = x;

    para = fcnParametersForSolver(p);
    
    
    T0 = para.T0;
    TR = para.TRw0*T0;
    mr = para.mr;
    Ti = para.Ti*T0;
    gamma = para.gamma;
    qr = para.qr/para.cv0;
    Rs = para.Rs; % Kushnir R*
    Us = para.Us;
    
    a1 = (gamma*Ti/TR + qr)/(gamma-Rs+qr);
    a2 = Us/(gamma-Rs+qr+1);
    b1= gamma-Rs+qr;
    
    rho = 1+mr;
    
    gx = (a1+a2.*rho + (T0/TR - a1-a2)*rho.^(-b1))*TR - 273.15 - Ttarget;
    
end 


function dgdx = dgdxfunction(gxfunction,x,gx) %#ok<DEFNU>
    % This function calculates numerically the jacobian dT/dx(2)
    % T: Maximal temperature during injection
    % x(1): Injection duration time
    % x(2): Injection rate
    
    eps = 0.001;
    delta = 1e-10;
    dx = x*eps+delta;
    xp = x + dx;
    gxp = gxfunction(xp);
    dgdx = (gxp-gx)/dx;
end
 


