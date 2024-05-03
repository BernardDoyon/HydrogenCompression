% ************************************************************************
%
%   main.m file¸
%
%           1- This main file solves the equations for the hydrogen gaz and 
%               rock temperatures during injection/storage/discharge. Two
%               methods of Injection/Discharge are possible
%
%               A) Imposing a function for dm/dt during storage/discharge 
%                           dm: variation of mass
%                           dt: variation of time 
%
%               B) Imposing a function for dP/dt during storage/discharge
%                           dP: variation of pressure
%                           dt: variation of time 
%           
%           2- The reservoir can be outside (surronded by air) or in the ground.
%               (vertical or horizontal). reservoir dimension and thermophysical
%               properties nedd to be specify in
%               fcnSelectPhysicalPropertiesOfCavity.m (edit this file !)
%
%           3- Heat transfer coefficient can be imposed or estimated from
%               Nusselt parameter. THIS IS NOT OK ... WE NEED TO IMPOSE A
%               MINMMAL VALUE FOR THE PRESSURE INSIDE THE INLET. IF WE USE
%               P INLET = P reservoir, WHEN P reservoir IS TO SMALL THIS PRODUCE 
%               INJECTION SPEED FOR H2 GAS TO HIGH (V of the order of 10^5 m/s !!!!)  
%
%           2- The numerical parameters for the simulations need to be
%               set from line 62
%
%           3- The injection/storage/scenario need to be specified 
%               The file InjectStoraExtractScenarioFunction.m can be
%               edited in order to wrtie the functions dm/dt or dP/dt. Some
%               default functions are coded for the moment (see line 119 of
%               this main file)
%
%           4- The tasks for the main file can be ajusted from line 164
%               (it's possible for instance to work on the figures without 
%                having to solve the equations if the solutions have already 
%                been calculated)   
%
%           5- The solvers are called from line 202
%
%           6- The code for the figures starts at line 279
%             
% ======================================================================
% Copyright (c) December 2023, Bernard Doyon (bdoyon@cegepgarneau.ca)
%
% ======================================================================       * 
%
%% ***********************************************************************

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
    % The file  fcnSelectPhysicalPropertiesOfCavity.m need to be edited
         % to ajust the size of the cavity and thermal properties of Wall
         % 3 exemples are given and can selected :
                 % 1- wall with 3 different layers
                 % 2- wall with 2 different layers (the second one is salt)
                 % 3- wall is salt 
                 % 4 Couteau Cavity
    [p.underground, p.geometry, p.Rw_int, p.H,p.Dinlet,p.ThermalPropertiesCavity] = fcnSelectPhysicalPropertiesOfCavity(1);  
    
    flagRealGas=1; %0- Simplified/Ideal gas  1- real gas model for Hydrogen
    
    %(if flagRealgas = 0 => p.model is used BUT p.model need to be set 
    p.model = 1; % 0 -> Ideal gaz: Z = 1
                 % 1 -> Simplified model: Z = Z0; dZdT = ZT0
    
    % The gas parameters (Hydrogen)
    p.R = 4.12448151675695e+03;         % Specific gas constant for hydrogen (J/Kg/K)
    p.Ti = 30 + 273.15;                 % Injected gas temperature (K)
    p.T0 = 30 + 273.15;                 % Initial gas temperature (K)
    p.P0 = 1e6;                         % Initial gas pressure (Pa)
    p.dz = p.H/2;                       % Vertical distance between inlet and middle of the reservoir (m)

    p.hc_int = 0;         % Heat transfer coefficient (H2 with interior wall)  (W/m2/K)
                           % if 0 => Heat transfer is estimated from
                           % Nusselt Number (see Couteau et al.
                           % Inertnational Journal of Hydrogen Energy 47
                           % (2022) 23060-23069
    
    if (p.underground == 0) % reservoir is outside; need to set air temperature and heat transfer coefficient between reservoir and air
       p.Tair = 30;
       p.hc_ext = 8;
    end
    p.TRw0 = p.T0;          % Initial rock temperature
  
    
    % Numerical solver options

    SolverOptions.solver = 1;       % 0- ode45 solver (Non-stiff explicit Runge Kutta)
                                    % 1- ode23s solver (stiff solver)

    SolverOptions.N = 50;           % Number of grid point in the reservoir
   
    SolverOptions.beta = 1.05;      % Beta parameter for the variable change r -> eta
    SolverOptions.Rp = 0;           % The minimum value of R in the rock where the temperature is not affected by heat conduction
                                    % When set to 0, eq A.7 of Kushnir 2012 is used
    SolverOptions.eps = 0.001;      % Parameter epsilon in eq A.7 of Kushnir
    SolverOptions.RelTol = 1e-4;    % Relative Tolerance for the numerical solver
    SolverOptions.AbsTol = 1e-6;    % Absolute Tolerance for the numerical solver
    SolverOptions.NtimeStep = 0;  % If NtimeStep = 0, the solution is return at every time step of the ode solver 
                                    % If NtimeStep = N, the solution is return for N time step of equal values.  
end
%
%% ********* END SET NUMERICAL PARAMETERS FOR SIMULATION ******************

%% ********************** SET INJECTION SCENARIO HERE *********************
%
for n= 1:1
 % The injection scenario parameters

    p.tp  = 12*3600;         % Period of the cycle  (s)
    p.mc = 2600/p.tp;        % Mass rate variation during injection/extraction (kg/s) (for dm/dt simulations)
   %p.mc = 0.06;
    p.pc = 118.5e6/p.tp;     % Pressure rate variation during injection/extraction (Pa/s) (for dP/dt simulations)

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
    Fie = 2;
    Gie = 2; 
end
%
%% ************************ END OF INJECTIOM SCENARIO *********************


%% ******** TASKS TO DO FOR THIS MAIN FILE ********************************
% 
% SOLVE_EQ -> numerical solutions are calculated 
% FIGURES - Plot options (to show the results) 
%               The figures are saves in "Figures" directory. 
%               This directory will be created in the parent directory
%
SOLVE_EQ = 1; %  %(0-> no; 1->yes)

% Figure 1 shows: the air temperature, the density and the pressure
% as a function of time for dm/dt and dP/dt imposed functions
% (Simple and complex gaz temperature are shown) 
FIGURE1 = 1; %(0->no; 1->yes)
Figure1Title = 'Comparing injection with constant dm/dt and constant dP/dt';
Figure1FileName = 'Figure1';

% Figure 2: Shows the temperature profile in the ground for 4 different
% time and for the case dm/dt = cte
FIGURE2 = 1; %(0-no; 1-> yes)
Figure2Title = {'Temperature profile in the ground as a function of time with dm/dt = cte',...
                '(Injection for 6 hours; Storage for 6 hours)'};
Figure2FileName = 'Figure2';         

FIGURE3 = 0;
Figure3Title = {'Temperature profile in the ground as a function of time with dm/dt = cte',...
                '(Injection for 6 hours; Storage for 6 hours)'};
Figure3FileName = 'Figure3';         

FIGURE4 = 0;
Figure4Title = {'Temperature profile in the ground as a function of time with dm/dt = cte',...
                '(Injection for 6 hours; Storage for 6 hours)'};
Figure4FileName = 'Figure4'; 

FIGURE5 = 1;

%% ************ END OF TAKS TO DO FOR THIS MAIN FILE **********************

%% ***************** SOLVE THE EQUATIONS  *********************************
%
if(SOLVE_EQ) % The numerical solutions are evaluated here
    for n =1:1
        
        if ~isfile('ScatterInterpFunctions.mat')
            fprintf('Generating interpolation functions and saving to ScatterInterpFunctions.mat \n');
            fcnGenerateHydrogenPropertiesFunctions()
        end
        
        [p.Fie, p.Gie] = fcnInjectStoraExtractScenarioFunctions(Fie,Gie); % function handlers for p.Fie and p.Gie
        
        p.Rw = p.ThermalPropertiesCavity(p.Rw_int,p.T0);
        p.kR = zeros(1,length(p.Rw));
        p.CpR = zeros(1,length(p.Rw));
        p.rhoR = zeros(1,length(p.Rw));
        kRT2 = zeros(1,length(p.Rw));
        % read thermal properties of wall and check if it depends on
        % temperature
        flagkR = 0;
        for i = 1:length(p.Rw)
            %thermal properties at initial temperature T0
            [p.Rw,  p.kR(i), p.CpR(i), p.rhoR(i)] = p.ThermalPropertiesCavity(p.Rw(i)-1e-6,p.T0);
            %thermal properties at 1.5*T0
            [p.Rw, kRT2(i), temp, temp] = p.ThermalPropertiesCavity(p.Rw(i)-1e-6,1.5*p.T0);
            if (p.kR(i) ~= kRT2(i))
               flagkR = 1; 
            end
        end    
       
        p.alpha = p.kR./p.CpR./p.rhoR; %alpha at initial temperature
        
        p.Ac_int = 2*pi*p.Rw_int*p.H;        % Surface of cylindrical cavity (m^2)
        p.V = pi*p.Rw_int^2*p.H;  
        
        para = fcnParametersForSolver(p); % Initialize the structure para
         
        % solving with function dm/dt: 
        tic; 
        fprintf('Solving for dm/dt function (Real gaz) \n');
        if (flagkR)
            [t1,rho1, Tair1, P1, Trock1,rrock1,Nc1] = fcnSolve_dmdt_kRvariable(NbCycles,para,SolverOptions,flagRealGas);
        else
            [t1,rho1, Tair1, P1, Trock1,rrock1,Nc1] = fcnSolve_dmdt(NbCycles,para,SolverOptions,flagRealGas); % Simulate with function dm/dt: 
        end
        time1= toc;
        
        fprintf('Solving for dP/dt function (Real gaz) \n');
        if (flagkR)
            [t2,rho2, Tair2, P2, Trock2,rrock2,Nc2] = fcnSolve_dPdt_kRvariable(NbCycles,para,SolverOptions,flagRealGas);    % Simulate with function dP/dt: 
        else
            [t2,rho2, Tair2, P2, Trock2,rrock2,Nc2] = fcnSolve_dPdt(NbCycles,para,SolverOptions,flagRealGas);    % Simulate with function dP/dt: 
        end
        time2 = toc;
        
        tscale = p.tp(end)/3600;
        t1plot = t1*tscale;
        t2plot = t2*tscale;
        
        Tair1plot = Tair1*para.T0-273.15;
        Tair2plot = Tair2*para.T0-273.15;
        
        rho0 = para.rho0;
        rho1plot = rho1*rho0;
        rho2plot = rho2*rho0;
     
        P1plot = P1*p.P0/1e6;
        P2plot = P2*p.P0/1e6;
       
        Trock1plot = Trock1*para.T0-273.15;
        Trock2plot = Trock2*para.T0-273.15;
    end
end
%
%% **************** END OF SOLVE THE EQUATIONS ****************************
 

%% ****************% CODE FOR FIGURES ******************************
%
for n=1:1
   
    currentFolder = pwd;
    FigurePath = 'Figures';
    if ~isfolder(FigurePath)
        mkdir(FigurePath)
    end
    %set size of figure on the screen 
    set(groot, 'defaultFigureUnits','normalized');
    set(groot, 'defaultFigurePosition',[0.2 0.2 0.8 0.8]);
    set(groot, 'defaultFigurePosition',[0 0 0.8 0.8]);

    if (FIGURE1)
        FigurePath = 'Figures';
        if ~isfolder(FigurePath)
            mkdir(FigurePath)
        end
        fig1= figure(1);
        subplot(3,1,1)
        plot(t1plot,Tair1plot,'--r',t2plot,Tair2plot,':b','Linewidth',2);
        legend('Nmercial model with dm/dt = cte','Numercial model with dP/dt = cte', 'Location','northeast');
        sgtitle(Figure1Title);
        xlabel('t (hours)');
        ylabel('Temperature (^oC)');
        ax = gca;
        ax.FontSize = 18; 
        % Plot solution for rho(t) for injection with dm/dt and dP/dt
        subplot(3,1,2)
        plot(t1plot,rho1plot,'--r',t2plot,rho2plot,':b','Linewidth',2);
        xlabel('t (hours)');
        ylabel('rho (kg/m^3)');
        ax = gca;
        ax.FontSize = 18;
        % Plot solution for P(t) for injection with dm/dt and dP/dt
        subplot(3,1,3)
        plot(t1plot,P1plot,'--r',t2plot,P2plot,':b','Linewidth',2);
        xlabel('t (hours)');
        ylabel('P (MPa)');
        ax = gca;
        ax.FontSize = 18;
        delete (strcat(currentFolder,'\',FigurePath,'\',Figure1FileName,'.png'));
        print(fig1,'-dpng', strcat(currentFolder,'\',FigurePath,'\',Figure1FileName));
    end
    
    if (FIGURE2)
        FigurePath = 'Figures'; %#ok<UNRCH>
        if ~isfolder(FigurePath)
            mkdir(FigurePath)
        end
       
        %set size of figure on the screen
        set(groot, 'defaultFigureUnits','normalized');
        set(groot, 'defaultFigurePosition',[0.2 0.2 0.8 0.8]);
        set(groot, 'defaultFigurePosition',[0 0 0.8 0.8]);
        
        fig2 = figure(2);
        
        snaptime = [.25 .5 .75 1];
        
        Tr1 = interp1(t1,Trock1plot,snaptime(1));
        plot(rrock1*p.Rw_int,Tr1,'--b','Linewidth',2);
        hold on
        Tr1 = interp1(t1,Trock1plot,snaptime(2));
        plot(rrock1*p.Rw_int,Tr1,':b','Linewidth',2);
        Tr1 = interp1(t1,Trock1plot,snaptime(3));
        plot(rrock1*p.Rw_int,Tr1,'-.b','Linewidth',2);
        Tr1 = interp1(t1,Trock1plot,snaptime(4));
        plot(rrock1*p.Rw_int,Tr1,'*','Linewidth',2);
        
     %   Tr1 = interp1(t2,Trock2plot,snaptime(4));
     %   plot(rrock2*p.Rw_int,Tr1,'-k','Linewidth',2);
        
        hold off
        
        for i = 1:(length(p.Rw)-1)
            %xline(p.Rw(i),'--r',{'First','casing'},'FontSize',24,'LineWidth',2);
            xline(p.Rw(i),'--r','FontSize',24,'LineWidth',2);
        end
        
        xlabel(strcat('Radius (m) (R_w=',sprintf(' %.2f',p.Rw_int),' m)'),'Fontsize',16);
        ylabel('Rock temperature (^oC)','Fontsize',16);
        
        txt1 = strcat('t =',int2str(snaptime(1)*tscale),' h');
        txt2 = strcat('t =',int2str(snaptime(2)*tscale),' h');
        txt3 = strcat('t =',int2str(snaptime(3)*tscale),' h');
        txt4 = strcat('t =',int2str(snaptime(4)*tscale),' h');
        
        legend(txt1,txt2,txt3,txt4,'Fontsize',14);
        
        xlim([p.Rw_int max(rrock1*p.Rw_int)]);
        
        title(Figure2Title,'Fontsize', 18);
  
        delete (strcat(currentFolder,'\',FigurePath,'\',Figure2FileName,'.png'));
        print(fig2,'-dpng', strcat(currentFolder,'\',FigurePath,'\',Figure2FileName));
    end
    
    if (FIGURE3)
        FigurePath = 'Figures'; %#ok<UNRCH>
        if ~isfolder(FigurePath)
            mkdir(FigurePath)
        end
        
        fig3 = figure(3);
        Tmin = ceil(min(min(Trock1plot)));
        Tmax = ceil(max(max(Trock1plot)));
        Rmin = p.Rw_int;
        Rmax = rrock1(end)*p.Rw_int;
        for t = 1:length(t1)
            plot(rrock1*p.Rw_int,Trock1plot(t,:),'LineWidth',2);
            ylim([Tmin Tmax]);
            xlim([Rmin Rmax]);
            xlabel(strcat('Radius (m) (R_w=',sprintf(' %.2f',p.Rw_int),' m)'));
            ylabel('Rock temperature (^oC)');
            title(Figure3Title,'FontSize',20);
            text(0.9,0.9,strcat('t = ',sprintf(' %.2f',t1plot(t)),'h'),'FontSize',24,'Units','normalized');
            ax = gca;
            ax.FontSize = 18;
              
            for i = 1:(length(p.Rw)-1)
                %xline(p.Rw(i),'--r',{'First','casing'},'FontSize',24,'LineWidth',2);
                xline(p.Rw(i),'--r','FontSize',24,'LineWidth',2);
            end
            M(t) = getframe(fig3);
        end
        v = VideoWriter(strcat(currentFolder,'\Figures\',Figure3FileName),'Uncompressed AVI'); %#ok<TNMLP>
        open(v)
        writeVideo(v,M)
        close(v)
    end
    
     if (FIGURE4)
        fig4 = figure(4); %#ok<UNRCH>
        Tmin = ceil(min(min(Trock1plot)));
        Tmax = ceil(max(max(Trock1plot)));
        v = linspace(Tmin, Tmax,100);
        y = [0 p.H];
        [X,Y] = meshgrid([0 0.99 rrock1]*p.Rw_int,y);
        for t = 1:length(t1)
            Z = [Tair1plot(t) Tair1plot(t) Trock1plot(t,:)];
            Z = [Z; Z];
            contourf(X,Y,Z,v,'LineStyle', 'none' );
            caxis([Tmin Tmax]);
            colormap(gray);
            c = colorbar;
            c.Label.String = 'Temperature (^oC)';
            c.Label.FontSize = 18;
            colormap(flipud(colormap));
            %xline(p.Rw.A,'--r',{'First','casing'},'FontSize',24,'LineWidth',5);
            %xline(p.Rw.B,'--b',{'Second','casing'},'FontSize',24,'LineWidth',5);
            set ( gca, 'ydir', 'reverse' );
            xlabel(strcat('Radius (m) (R_w=',sprintf(' %.2f',p.Rw_int),' m)'));
            ylabel('Depth of the well (m)')
           % title(Figure3Title,'FontSize',20);
            text(0.9,0.9,strcat('t = ',sprintf(' %.1f',t1plot(t)),'h'),'FontSize',24,'Units','normalized');
            if (para.Fie(t)>0)
                text(0.9,0.85,'(Injection)','FontSize',24,'Units','normalized');
            elseif (para.Fie(t)<0)
                text(0.9,0.85,'(Discharge)','FontSize',24,'Units','normalized');
            else
                text(0.9,0.85,'(Storage)','FontSize',24,'Units','normalized');
            end
            ax = gca;
            ax.FontSize = 18;
             for i = 1:(length(p.Rw)-1)
                %xline(p.Rw(i),'--r',{'First','casing'},'FontSize',24,'LineWidth',2);
                xline(p.Rw(i),'--r','FontSize',24,'LineWidth',2);
            end
            M(t) = getframe(fig4);
        end
         v = VideoWriter(strcat(currentFolder,'\Figures\',Figure4FileName),'Uncompressed AVI'); %#ok<TNMLP>
         open(v)
         writeVideo(v,M)
         close(v)
    end
    
    
    if (FIGURE5) % 
        fig5 = figure(5);
        t = t1;
        T = Tair1*p.T0;
        P = P1plot*1e6;
        rho = rho1plot;
        TRw = Trock1(:,1)*para.T0;

       
        Z = para.Z_TP(T,P); 
        Ti = ones(length(t),1)*para.Ti*p.T0;
     
        mdot2 = diff(rho)*para.V./(diff(t)*p.tp);
        mdot = [mdot2(1) mdot2']';
        
        %density of injected gaz;
        rho_i = Z.*rho.*T./para.Z_TP(Ti,P)./Ti;
        % velocity of injected gaz
        
        v = mdot./rho_i./pi/(para.Dinlet/2)^2;
   
        mu_inlet = para.mu(Ti,P);
        mu = para.mu(T,P);
        kH2 = para.kH2(T,P);
        cp = para.Cp(T,P);
        Beta = 1./T + 1./Z.*para.dZdT_P(T,P); 
        
        Re_in = 4*mdot./mu_inlet/pi/para.Dinlet;
        Ra = 9.8*Beta.*cp.*rho.^2*para.L^3.*(T-TRw)./mu./kH2;
    
        a = 0.32;
        Nu_for = a*Re_in.^0.67;
        Nu_nat = 0.104*Ra.^0.352;
        Ratio = Nu_for./Nu_nat;
        
        if (p.hc_int==0)
            hc = kH2.*(Nu_for + Nu_nat)/para.L;
        else
            hc = ones(length(t),1)*p.hc_int;
        end
        
        Q = -hc.*p.Ac_int.*(T-TRw); %kW
        
        hinlet = para.H(Ti,P);
        hgaz = para.H(T,P);
        Term3 = mdot.*(hinlet-hgaz+v.^2/2);
        Term1 = para.V*Beta.*T.*Gie.*p.pc;
        
        subplot(2,1,1)
        plot(t1plot,v);
        xlabel('t (hours)');
        ylabel('Injected gas speed (m/s))');
        subplot(2,1,2);
        plot(t1plot,hc);
        xlabel('t (hours)');
        ylabel('Heat transfer coefficient (W/m^2/s)');
   
        
        
       
     %   delete (strcat(currentFolder,'\',FigurePath,'\',Figure4FileName,'.png'));
     %   print (fig4,'-dpng', strcat(currentFolder,'\',FigurePath,'\',Figure4FileName));  
    end
    
    
    
    
   
    
    
    set(groot, 'defaultFigureUnits','default');
    set(groot, 'defaultFigurePosition','default');
end % The code to plot the figures
%
%% ************** END OF CODE FOR FIGURES ***************************
