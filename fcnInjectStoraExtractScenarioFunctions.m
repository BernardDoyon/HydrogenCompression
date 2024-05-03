function [FunctFie,FunctGie] = fcnInjectStoraExtractScenarioFunctions(Fie,Gie)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This files select the hander for de function Fie(t) and and Gie(t): 
%
%   If we inject/discharge by imposing a function for the mass variation 
%   then we have:
%           
%           dm/dt = m_c*F_ie(t) 
%                                       
%           with m_c = constant (kg/s) (to be set in the main file with the strucutre p.mc)
%                F_ie(t) a function of time (positive for injection, negative for discharge, 0 for storage)
%
%   If we injecte/discharge by imposing a function for the pressur variation 
%   then we have:
%           
%           dp/dt = p_c*G_ie(t) 
%                                       
%           with p_c = constant (Pa/s) (to be set in the main file with the strucutre p.pc)
%                G_ie(t) a function of time (positive for injection, negative for discharge, 0 for storage)
%
%           Remark: We set Gie =0 for storage only to let the algorithm know that there is a storage phase. 
%                    Of course the pressure will change during storage
%                    (because the change in temperature)
%
%    INPUT:
%           Fie: Integer to select the function F_ie(t)
%                   1 - Constant injection rate from t*= 0 to t* = 1
%                   2 - Constant injection rate from t* = 0 to t* = 0.5 and
%                       storage from, t* = 0.5 to t* = 1
%                   3 - Constant injection rate from t* = 0 to t* = 0.25; 
%                       storage from t* = 0.25 to t* = 0.5;
%                       discharge from t* = 0.5 to t* = 0.75;
%                       storage frome t* = 0.75 to t* = 1;
%                   4 - A specific test with a decrease of the injection rate
%
%           Gie: Integer to select the function G_ie(t)
%                   1 - Constant positive pressure rate from t*= 0 to t* = 1
%                   2 - Constant positive pressure rate from t* = 0 to t* = 0.5 and
%                       storage fro, t* = 0.5 to t* = 1
%                   3 - Constant positive pressure rate from t* = 0 to t* = 0.25; 
%                       storage from t* = 0.25 to t* = 0.5;
%                       constant negative pressure rate from t* = 0.5 to t* = 0.75;
%                       storage frome t* = 0.75 to t* = 1;
%                   4 - Couteau scenario 
%                           i) Constant dP/dt for 120 s and
%                           ii) dP/dt=0 from 120s to 140 s
%
%   OUTPUT
%
%           FunctFie: Handler to function F_ie(t)
%
%           FunctGie: Handler to function G_ie(t)
%


    switch(Fie)
        case 1
             FunctFie = @funcFieInjectionOnly;
        case 2
             FunctFie = @funcFieInjectionStorage;
        case 3
             FunctFie = @funcFieCompleteCycle;
        case 4
             FunctFie = @funcFieSpecificTest;
        case 5
             FunctFie = @funcFieStorageOnly;
    end
    
    switch(Gie) 
        case 1
             FunctGie = @funcGieInjectionOnly;
        case 2
             FunctGie = @funcGieInjectionStorage;
        case 3
             FunctGie = @funcGieCompleteCycle;
        case 4
             FunctGie = @funcGieCouteauScenario;
        case 5
            temp = readtable('CouteauFigure3PressionInlet.xlsx');
            temp = table2array(temp);
            tP = temp(:,1);
            P =temp(:,2); % read Couteau inlet Pressure
            tPs =tP(1:26); %
            Ps = P(1:26); %

            pp = polyfit(tPs,Ps,2); % smooth the data for injection (t< 131 s)
            tfit = linspace(0,131.3,100);
            Pfit = polyval(pp,tfit);
            
            dPdt = diff(Pfit)./diff(tfit); % estimate dP/dt for injection
            tdPdt = (tfit(2:end)+tfit(1:(end-1)))/2;
            meandPdt = mean(dPdt);
            dPdt = [dPdt dPdt(end)];
            tdPdt = [tdPdt 131.3];
            dPdt = [dPdt zeros(1,9)]/meandPdt;
            tdPdt = [tdPdt 132 133 134 135 136 137 138 139 140]/140;
            FunctGie = griddedInterpolant(tdPdt,dPdt);
        case 6
           FunctGie = @funcGieCouteauScenarioFigure6; 
    end
end

function Fie = funcFieSpecificTest(t)
    
        mr = 314.4041;
        t1 = 1/4;
        t2 = 1/2;
        t3 = 3/4;
     
        % the extracted mass
        if (t <= t1)
            Fie = 1;
        elseif (t <= t2)
            Fie=0;
        elseif (t <= t3)
            Fie = (3*766.5*t^2 - 2*1652.9*t + 863.5)/mr;
        else
            Fie=0;
        end
        
end


function Fie = funcFieInjectionOnly(t)
        Fie = 1;
end


function Fie = funcFieInjectionStorage(t)
    if (t<=0.5)
        Fie = 1;
    else
        Fie=0;
    end
end

function Fie = funcFieCompleteCycle(t)
    t1 = 1/4;
    t2 = 1/2;
    t3 = 3/4;
    CD = t1/(t3-t2); 
        % For this exemple, we impose the injected mass to be equal to
        % the extracted mass
        % It is not necessary BUT the injected mass has to be greater then
        % the extracted mass 
    if (t <= t1)
        Fie = 1;
    elseif (t <= t2)
        Fie=0;
    elseif (t <= t3)
        Fie= -CD;
    else
        Fie=0;
    end
end


function Fie = funcFieStorageOnly(t)
    Fie = 0;
end

function Gie = funcGieInjectionOnly(t) %#ok<*INUSD>
        Gie = 1;
end


function Gie = funcGieInjectionStorage(t)
    if (t<=0.5)
        Gie = 1;
    else
        Gie=0;
    end
end

function Gie = funcGieCompleteCycle(t)
    t1 = 1/4;
    t2 = 1/2;
    t3 = 3/4;
    if (t <= t1)
        Gie = 1;
    elseif (t <= t2)
        Gie = 0;
    elseif (t <= t3)
        Gie = -0.9;
    else
        Gie = 0;
    end
end


function Gie = funcGieCouteauScenario(t)
    t1 = 120/140;
    if (t <= t1)
        Gie = 1;
    else
        Gie=0;
    end
end


function Gie = funcGieCouteauScenarioFigure6(t)
    t1 = 120/140;
    if (t <= t1)
        Gie = 1;
    else
        Gie=0;
    end
end
