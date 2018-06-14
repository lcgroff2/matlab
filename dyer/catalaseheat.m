% Calculate Catalase heme B heating after disproportionation of H2O2 
% assuming that all of the heat is transfered into the product molecule:

T1 = 298;  %(K) - Assume room temp to start.

delta_H = 100; %(kJ/mol) - Estimated total heat released during catalytic 
               %turnover from Bustamante, et. al.

Cp_ZnTPP = 1.35/1000; %(kJ g-1 K-1) - per gram heat capacity of ZnTPP
MW_ZnTPP = 678.11;   %(g/mol) - molecular weight of ZnTPP
Cpmol_ZnTPP = Cp_ZnTPP*MW_ZnTPP; %(kJ mol-1 K-1) - molar heat capacity of 
                                 %ZnTPP

dTZn = delta_H/Cpmol_ZnTPP; %(K) - dH = CpdT rearranged to isolate dT, using 
                     %remaining heat to calculate temperature change of 
                     %water vapor.
T_final = T1+dTZn; %(K) - Add up temperature changes from total heat 
                   %transfer.