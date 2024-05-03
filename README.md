# HydrogenCompression

(See MatlabComputerProgramSummary.docx to visualize the dependencies of the different functions)
 
 main.m file¸

           1- This main file solves the equations for the hydrogen gaz and 
               rock temperatures during injection/storage/discharge. Two
               methods of Injection/Discharge are possible

               A) Imposing a function for dm/dt during storage/discharge 
                           dm: variation of mass
                           dt: variation of time 

               B) Imposing a function for dP/dt during storage/discharge
                           dP: variation of pressure
                           dt: variation of time 
           
           2- The reservoir can be outside (surronded by air) or in the ground.
               (vertical or horizontal). reservoir dimension and thermophysical
               properties nedd to be specify in
               fcnSelectPhysicalPropertiesOfCavity.m (edit this file !)

           3- Heat transfer coefficient can be imposed or estimated from
               Nusselt parameter.  

           2- The numerical parameters for the simulations need to be
               set from line 62

           3- The injection/storage/scenario need to be specified 
               The file InjectStoraExtractScenarioFunction.m can be
               edited in order to wrtie the functions dm/dt or dP/dt. Some
               default functions are coded for the moment (see line 119 of
               this main file)

           4- The tasks for the main file can be ajusted from line 164
               (it's possible for instance to work on the figures without 
                having to solve the equations if the solutions have already 
                been calculated)   

           5- The solvers are called from line 202

           6- The code for the figures starts at line 279
             
 ======================================================================
 Copyright (c) December 2023, Bernard Doyon (bdoyon@cegepgarneau.ca)
 ======================================================================     
