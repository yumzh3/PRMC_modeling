# modeling for periodically replenished magma chambers, fractionaly crystallization is used in modeling
# calculate steady state liquid compositions after crystallization and after mixing for two type of sequences.
# 1) crystallization-eruption-mixing (fem), with a given (1-F)/E ratio with different (1-F) values
# steady state liquids after crystallization are saved in the data frame "steady_afterxtal_PRMC_fem_combined"
# corresponding liquids after mixing are saved in the data frame "steady_aftermix_PRMC_fem_combined"
# liquids path to each steady state liquid after crystallization are saved in the dictionary "ALL_DATA_afterxtal_PRMC"
# liquids path to corresponding liquid after mixing for each steady state are saved in the dictionary "ALL_DATA_aftermix_PRMC"
# crystallization proportion, eruption proportion and the refill proportion are saved in the data frame "FEM_df".
# 2) crystallization-mixing-eruption (fme), with a given (1-F)/E ratio with different (1-F) values
# steady state liquids after mixing are saved in the data frame "steady_aftermix_PRMC_fme_combined"
# corresponding liquids after crystallization are saved in the data frame "steady_afterxtal_PRMC_fme_combined"
# liquids path to each steady state liquid after mixing are saved in the dictionary "ALL_DATA_aftermix_PRMC"
# liquids path to corresponding liquid after crystallization for each steady state are saved in the dictionary "ALL_DATA_afterxtal_PRMC"
# crystallization proportion, eruption proportion and the refill proportion are saved in the data frame "FME_df".
# also modeling for pure fractional and equilibrium crystallization
# calculate the liquid fraction, phase fractions, liquid and phase compositions in wt.% per each step during crystallization
# calculte either fractional or equilibrium crystallization
# results are saved in data frames "LLD_df_frac" and "LLD_df_equ"
# last modified: Sep 26, 2025

from wl1989stoich_2023 import *
from wl1989kdcalc_2023 import *
from wl1989models_2023 import *
from wlState_2023 import *
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import random


# %% load starting compositions
LLDinput = pd.read_csv('local location/parental_magma_data.csv')
pm_location = 'OJ_MORBMgO10_final'  # can add more lines of different starting compositions in the data file
index = np.where(LLDinput['location'] == pm_location)[0][0]
P = 1 # pressure for crystallization, unit is bar  

# %% PRMC modeling
fxtale_ratio = 3  # (1-F)/E ratio, crystallization/eruption ratio for each cycle
# F_xtal_PRMC = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75]  # (1-F), crystallization fraction for each cycle, cannot be higher than 1
F_xtal_PRMC = [0.05]
F_erupt_PRMC = [x/fxtale_ratio for x in F_xtal_PRMC]  # E, eruption fraction for each cycle, constrained by (1-F)/E with a given (1-F)
ALL_DATA_afterxtal_PRMC = {}
ALL_DATA_aftermix_PRMC = {}

## == PRMC: crystallization-eruption-mixing (fem)  
for simulation_number in range(0,len(F_xtal_PRMC)):     
    parental_magma = LLDinput.loc[index,['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']]  # starting composition for the first cycle and is one of the endmember for mixing in each cycle
    parental_magma['FeO'] = parental_magma['FeO']*0.9  # assume FeO/FeO(t) = 0.9 for MORB
    parental_magma_trace = LLDinput.loc[index,['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']]
    f_xtal_PRMC = F_xtal_PRMC[simulation_number]
    f_erupt_PRMC = F_erupt_PRMC[simulation_number]
    f_replenish_PRMC = f_xtal_PRMC + f_erupt_PRMC  # (1-F)+E, refill fraction for each cycle
    f_target = 1-f_xtal_PRMC  # targetd liquid fraction in each crystallization cycle
    magma = LLDinput.loc[index,['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']]
    magma['FeO'] = magma['FeO']*0.9  # assume FeO/FeO(t) = 0.9 for MORB
    magma_trace = LLDinput.loc[index,['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']]
    magma_afterxtal_PRMC_fem = {element: [] for element in ['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO',\
                                                            'La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']}
    magma_PRMC_fem = {element: [] for element in ['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO',\
                                                  'La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']}
    ## determine the iteration number based on refill fraction
    if f_replenish_PRMC >= 0.75:
        PRMC_iter_max = 30
    elif f_replenish_PRMC >= 0.5 and f_replenish_PRMC < 0.75:
        PRMC_iter_max = 60
    elif f_replenish_PRMC >= 0.1 and f_replenish_PRMC < 0.5:
        PRMC_iter_max = 200 
    elif f_replenish_PRMC >= 0.05 and f_replenish_PRMC < 0.1:
        PRMC_iter_max = 300
    elif f_replenish_PRMC >= 0.03 and f_replenish_PRMC < 0.05:
        PRMC_iter_max = 400
    elif f_replenish_PRMC >= 0.01 and f_replenish_PRMC < 0.03:
        PRMC_iter_max = 700
    elif f_replenish_PRMC >= 0.003 and f_replenish_PRMC < 0.01:
        PRMC_iter_max = 1000
    elif f_replenish_PRMC >= 0.001 and f_replenish_PRMC < 0.003:
        PRMC_iter_max = 1500
    elif f_replenish_PRMC < 0.001:
        PRMC_iter_max = 2000        
    print('fem_iter_max is '+str(PRMC_iter_max))
    
    for PRMC_iter in range(0,PRMC_iter_max):  # start the iteration
        ## crystallization using fractional crystallization
        system_components = magma.to_dict()
        T_system_components = oxideToComponent(system_components)
        t_start = get_first_T(T_system_components, P, kdCalc = kdCalc_langmuir1992)
        tstep = 0.5
        trace_components = magma_trace.to_dict()
        fl,fa_dict,major_oxide_dict,major_phase_oxide_dict,trace_dict,T_LLD  = frac_model_trange_PRMC(f_target, t_start, tstep, system_components,trace_components,kd_trace,P,kdCalc = kdCalc_langmuir1992)
        T_LLD_dict = {'T_C':T_LLD}
        T_df = pd.DataFrame(T_LLD_dict)
        T_df = T_df-273.15  # temperature in Celsius degree
        fl_dict = {'fl':fl}
        fl_df = pd.DataFrame(fl_dict)  # liquid fraction during crystallization
        fa_df = pd.DataFrame(fa_dict)  # mineral phase fractions during crystallization
        major_oxide_df = pd.DataFrame(major_oxide_dict)
        major_ol_oxide_df = pd.DataFrame(major_phase_oxide_dict['ol'])
        major_cpx_oxide_df = pd.DataFrame(major_phase_oxide_dict['cpx'])
        major_plg_oxide_df = pd.DataFrame(major_phase_oxide_dict['plg'])
        trace_dict_df = pd.DataFrame(trace_dict)
        ## liquid and mineral compositions for crystallization 
        LLD_df = pd.concat([T_df,fl_df,fa_df,major_oxide_df,major_ol_oxide_df,major_cpx_oxide_df,major_plg_oxide_df,trace_dict_df],axis=1)
        LLD_df.columns = ['T_C','f_liq','f_plg','f_cpx','f_ol','liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
                        'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','olSiO2',\
                            'olTiO2','olAl2O3','olFeO','olMgO','olK2O','olMnO','olNa2O','olP2O5','olCaO','olNiO',\
                                'cpxSiO2','cpxTiO2','cpxAl2O3','cpxFeO','cpxMgO','cpxK2O','cpxMnO','cpxNa2O',\
                                    'cpxP2O5','cpxCaO','cpxNiO','plgSiO2','plgTiO2','plgAl2O3','plgFeO','plgMgO',\
                                        'plgK2O','plgMnO','plgNa2O','plgP2O5','plgCaO','plgNiO','liq_La','liq_Ce',\
                                            'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                                                'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
        LLD_df['liq_FeOt'] = LLD_df['liq_FeO']/0.9  # assume FeO/FeO(t) = 0.9 for MORB
        LLD_df['Fo'] = 100/(1+LLD_df['olFeO']/LLD_df['olMgO']*40.3/71.84)
        LLD_df['olNippm'] = LLD_df['olNiO']*58.6934/74.69*10**4
        LLD_df['olMnppm'] = LLD_df['olMnO']*54.938/70.94*10**4
        LLD_df['liq_Ni'] = LLD_df['liq_NiO']*58.6934/74.69*10**4
        LLD_df['liq_FeOtMnO'] = LLD_df['liq_FeOt']/LLD_df['liq_MnO']
        
        f_var = list(abs(LLD_df['f_liq']-(1-f_xtal_PRMC)))
        liq_index_PRMC = f_var.index(min(f_var))  # get the liquid that reaches to the targeted crystallization proportion, this is the other endmember in the mixing
        liq_afterxtal = LLD_df.loc[liq_index_PRMC,['liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
                        'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO']]
        liq_afterxtal.index = ['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']
        ## mixing, starting compositions for the next PRMC iteration
        ## magma = [1-(1-F+E)]*liquid after crystallization + (1-F+E)*parental magma
        magma = liq_afterxtal*(1-f_xtal_PRMC-f_erupt_PRMC)+parental_magma*f_replenish_PRMC  # major components
        liq_trace_afterxtal = LLD_df.loc[liq_index_PRMC,['liq_La','liq_Ce',\
            'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']]
        liq_trace_afterxtal.index = ['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']
        magma_trace = liq_trace_afterxtal*(1-f_xtal_PRMC-f_erupt_PRMC)+parental_magma_trace*f_replenish_PRMC  # trace components
        for element in ['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']:
            magma_afterxtal_PRMC_fem[element].append(liq_afterxtal[element])  # liquid path to the steady state after crystallization
            magma_PRMC_fem[element].append(magma[element])  # liquid path to the corresponding liquids after mixing for the steady state
        for element in ['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']:
            magma_afterxtal_PRMC_fem[element].append(liq_trace_afterxtal[element])
            magma_PRMC_fem[element].append(magma_trace[element])
        print(PRMC_iter)
    ## create data frames for saving results        
    magma_afterxtal_PRMC_fem = pd.DataFrame(magma_afterxtal_PRMC_fem)
    magma_PRMC_fem = pd.DataFrame(magma_PRMC_fem)
    magma_afterxtal_PRMC_fem.columns = ['liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
                    'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','liq_La','liq_Ce',\
                        'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                            'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
    magma_PRMC_fem.columns = ['liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
                    'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','liq_La','liq_Ce',\
                        'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                            'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
    magma_afterxtal_PRMC_fem['liq_FeOt'] = magma_afterxtal_PRMC_fem['liq_FeO']/0.9
    magma_PRMC_fem['liq_FeOt'] = magma_PRMC_fem['liq_FeO']/0.9
    magma_afterxtal_PRMC_fem['liq_Ni'] = magma_afterxtal_PRMC_fem['liq_NiO']*58.6934/74.69*10**4
    magma_PRMC_fem['liq_Ni'] = magma_PRMC_fem['liq_NiO']*58.6934/74.69*10**4
    magma_erupted_fem = magma_afterxtal_PRMC_fem
    key_name_afterxtal_PRMC_fem = str(simulation_number)+'_fem_afterxtal'
    key_name_aftermix_PRMC_fem = str(simulation_number)+'_fem_aftermix'
    ALL_DATA_afterxtal_PRMC[key_name_afterxtal_PRMC_fem] = magma_afterxtal_PRMC_fem  # liquid paths to the steady state after crystallization for different (1-F) values
    ALL_DATA_aftermix_PRMC[key_name_aftermix_PRMC_fem] = magma_PRMC_fem  # liquid paths to the corresponding liquids after mixing for the steady state for different (1-F) values
    print('simulation number is '+str(simulation_number))  
## collect and compiled the steady state liquid compositions after crystallization and corresponding liquid compositions after mixing
steady_aftermix_PRMC = []
for simulation_number in range(0,len(F_xtal_PRMC)):
    key_name_aftermix_PRMC = str(simulation_number)+'_fem_aftermix'
    steady_aftermix_PRMC.append(ALL_DATA_aftermix_PRMC[key_name_aftermix_PRMC].tail(1))
steady_aftermix_PRMC_fem_combined = pd.concat(steady_aftermix_PRMC)
steady_afterxtal_PRMC = []
for simulation_number in range(0,len(F_xtal_PRMC)):
    key_name_afterxtal_PRMC = str(simulation_number)+'_fem_afterxtal'
    steady_afterxtal_PRMC.append(ALL_DATA_afterxtal_PRMC[key_name_afterxtal_PRMC].tail(1))
steady_afterxtal_PRMC_fem_combined = pd.concat(steady_afterxtal_PRMC)
FEM_refill = [x+y for x,y in zip(F_xtal_PRMC,F_erupt_PRMC)]
FxtalE_ratio = [fxtale_ratio]*len(F_xtal_PRMC)
FEM_df = pd.DataFrame({'xtal=1-F':F_xtal_PRMC,'erupt=E':F_erupt_PRMC,'refill=1-F+E':FEM_refill,'(1-F)/E':FxtalE_ratio})    

## == PRMC: crystallization-mixing-eruption (fme)  
for simulation_number in range(0,len(F_xtal_PRMC)):       
    parental_magma = LLDinput.loc[index,['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']]  # starting composition for the first cycle and is one of the endmember for mixing in each cycle
    parental_magma['FeO'] = parental_magma['FeO']*0.9  # assume FeO/FeO(t) = 0.9 for MORB
    parental_magma_trace = LLDinput.loc[index,['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']]
    f_xtal_PRMC = F_xtal_PRMC[simulation_number]
    f_erupt_PRMC = F_erupt_PRMC[simulation_number]
    f_replenish_PRMC = f_xtal_PRMC + f_erupt_PRMC  # (1-F)+E, refill fraction for each cycle
    f_target = 1-f_xtal_PRMC
    magma = LLDinput.loc[index,['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']]
    magma['FeO'] = magma['FeO']*0.9  # assume FeO/FeO(t) = 0.9 for MORB
    magma_trace = LLDinput.loc[index,['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']]
    magma_afterxtal_PRMC_fme = {element: [] for element in ['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO',\
                                                            'La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']}
    magma_PRMC_fme = {element: [] for element in ['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO',\
                                                  'La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']}
    ## determine the iteration number based on refill fraction
    if f_replenish_PRMC >= 0.75:
        PRMC_iter_max = 30
    elif f_replenish_PRMC >= 0.5 and f_replenish_PRMC < 0.75:
        PRMC_iter_max = 60
    elif f_replenish_PRMC >= 0.1 and f_replenish_PRMC < 0.5:
        PRMC_iter_max = 200 
    elif f_replenish_PRMC >= 0.05 and f_replenish_PRMC < 0.1:
        PRMC_iter_max = 300
    elif f_replenish_PRMC >= 0.03 and f_replenish_PRMC < 0.05:
        PRMC_iter_max = 400
    elif f_replenish_PRMC >= 0.01 and f_replenish_PRMC < 0.03:
        PRMC_iter_max = 700
    elif f_replenish_PRMC >= 0.003 and f_replenish_PRMC < 0.01:
        PRMC_iter_max = 1000
    elif f_replenish_PRMC >= 0.001 and f_replenish_PRMC < 0.003:
        PRMC_iter_max = 1500
    elif f_replenish_PRMC < 0.001:
        PRMC_iter_max = 2000        
    print('fme_iter_max is '+str(PRMC_iter_max))
    
    for PRMC_iter in range(0,PRMC_iter_max):  # start the iteration
        ## crystallization using fractional crystallization
        system_components = magma.to_dict()
        T_system_components = oxideToComponent(system_components)
        t_start = get_first_T(T_system_components, P, kdCalc = kdCalc_langmuir1992)
        tstep = 0.5
        trace_components = magma_trace.to_dict()
        fl,fa_dict,major_oxide_dict,major_phase_oxide_dict,trace_dict,T_LLD  = frac_model_trange_PRMC(f_target, t_start, tstep, system_components,trace_components,kd_trace,P,kdCalc = kdCalc_langmuir1992)
        T_LLD_dict = {'T_C':T_LLD}
        T_df = pd.DataFrame(T_LLD_dict)
        T_df = T_df-273.15  # temperature in Celsius degree
        fl_dict = {'fl':fl}
        fl_df = pd.DataFrame(fl_dict)  # liquid fraction during crystallization
        fa_df = pd.DataFrame(fa_dict)  # mineral phase fractions during crystallization
        major_oxide_df = pd.DataFrame(major_oxide_dict)
        major_ol_oxide_df = pd.DataFrame(major_phase_oxide_dict['ol'])
        major_cpx_oxide_df = pd.DataFrame(major_phase_oxide_dict['cpx'])
        major_plg_oxide_df = pd.DataFrame(major_phase_oxide_dict['plg'])
        trace_dict_df = pd.DataFrame(trace_dict)
        ## liquid and mineral compositions for crystallization 
        LLD_df = pd.concat([T_df,fl_df,fa_df,major_oxide_df,major_ol_oxide_df,major_cpx_oxide_df,major_plg_oxide_df,trace_dict_df],axis=1)
        LLD_df.columns = ['T_C','f_liq','f_plg','f_cpx','f_ol','liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
                        'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','olSiO2',\
                            'olTiO2','olAl2O3','olFeO','olMgO','olK2O','olMnO','olNa2O','olP2O5','olCaO','olNiO',\
                                'cpxSiO2','cpxTiO2','cpxAl2O3','cpxFeO','cpxMgO','cpxK2O','cpxMnO','cpxNa2O',\
                                    'cpxP2O5','cpxCaO','cpxNiO','plgSiO2','plgTiO2','plgAl2O3','plgFeO','plgMgO',\
                                        'plgK2O','plgMnO','plgNa2O','plgP2O5','plgCaO','plgNiO','liq_La','liq_Ce',\
                                            'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                                                'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
        LLD_df['liq_FeOt'] = LLD_df['liq_FeO']/0.9  # assume FeO/FeO(t) = 0.9 for MORB
        LLD_df['Fo'] = 100/(1+LLD_df['olFeO']/LLD_df['olMgO']*40.3/71.84)
        LLD_df['olNippm'] = LLD_df['olNiO']*58.6934/74.69*10**4
        LLD_df['olMnppm'] = LLD_df['olMnO']*54.938/70.94*10**4
        LLD_df['liq_Ni'] = LLD_df['liq_NiO']*58.6934/74.69*10**4
        LLD_df['liq_FeOtMnO'] = LLD_df['liq_FeOt']/LLD_df['liq_MnO']
        
        f_var = list(abs(LLD_df['f_liq']-(1-f_xtal_PRMC)))
        liq_index_PRMC = f_var.index(min(f_var))  # get the liquid that reaches to the targeted crystallization proportion, this is the other endmember in the mixing
        liq_afterxtal = LLD_df.loc[liq_index_PRMC,['liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
                        'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO']]
        liq_afterxtal.index = ['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']
        ## mixing, starting compositions for the next PRMC iteration
        ## magma = {[1-(1-F)]*liquid after crystallization + (1-F+E)*parental magma}/[1-(1-F)+(1-F+E)]
        magma = (liq_afterxtal*(1-f_xtal_PRMC)+parental_magma*f_replenish_PRMC)/(1-f_xtal_PRMC+f_replenish_PRMC)  # major components
        liq_trace_afterxtal = LLD_df.loc[liq_index_PRMC,['liq_La','liq_Ce',\
            'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']]
        liq_trace_afterxtal.index = ['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']
        magma_trace = (liq_trace_afterxtal*(1-f_xtal_PRMC)+parental_magma_trace*f_replenish_PRMC)/(1-f_xtal_PRMC+f_replenish_PRMC)  # trace components
        for element in ['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']:
            magma_afterxtal_PRMC_fme[element].append(liq_afterxtal[element])  # liquid path to the corresponding liquids after crystallization for the steady state
            magma_PRMC_fme[element].append(magma[element])  # liquid path to the steady state liquids after mixing
        for element in ['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']:
            magma_afterxtal_PRMC_fme[element].append(liq_trace_afterxtal[element])
            magma_PRMC_fme[element].append(magma_trace[element])
        print(PRMC_iter)
    ## create data frames for saving results        
    magma_afterxtal_PRMC_fme = pd.DataFrame(magma_afterxtal_PRMC_fme)
    magma_PRMC_fme = pd.DataFrame(magma_PRMC_fme)
    magma_afterxtal_PRMC_fme.columns = ['liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
                    'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','liq_La','liq_Ce',\
                        'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                            'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
    magma_PRMC_fme.columns = ['liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
                    'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','liq_La','liq_Ce',\
                        'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                            'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
    magma_afterxtal_PRMC_fme['liq_FeOt'] = magma_afterxtal_PRMC_fme['liq_FeO']/0.9
    magma_PRMC_fme['liq_FeOt'] = magma_PRMC_fme['liq_FeO']/0.9
    magma_afterxtal_PRMC_fme['liq_Ni'] = magma_afterxtal_PRMC_fme['liq_NiO']*58.6934/74.69*10**4
    magma_PRMC_fme['liq_Ni'] = magma_PRMC_fme['liq_NiO']*58.6934/74.69*10**4
    magma_erupted_fme = magma_PRMC_fme
    key_name_afterxtal_PRMC_fme = str(simulation_number)+'_fme_afterxtal'
    key_name_aftermix_PRMC_fme = str(simulation_number)+'_fme_aftermix'
    ALL_DATA_afterxtal_PRMC[key_name_afterxtal_PRMC_fme] = magma_afterxtal_PRMC_fme  # liquid paths to the corresponding liquids after crystallization for the steady state for different (1-F) values
    ALL_DATA_aftermix_PRMC[key_name_aftermix_PRMC_fme] = magma_PRMC_fme  # liquid paths to the steady state after mixing for different (1-F) values
    print('simulation number is '+str(simulation_number)) 
## collect and compiled the steady state liquid compositions after mixing and corresponding liquid compositions after crystallization    
steady_aftermix_PRMC = []
for simulation_number in range(0,len(F_xtal_PRMC)):
    key_name_aftermix_PRMC = str(simulation_number)+'_fme_aftermix'
    steady_aftermix_PRMC.append(ALL_DATA_aftermix_PRMC[key_name_aftermix_PRMC].tail(1))
steady_aftermix_PRMC_fme_combined = pd.concat(steady_aftermix_PRMC)
steady_afterxtal_PRMC = []
for simulation_number in range(0,len(F_xtal_PRMC)):
    key_name_afterxtal_PRMC = str(simulation_number)+'_fme_afterxtal'
    steady_afterxtal_PRMC.append(ALL_DATA_afterxtal_PRMC[key_name_afterxtal_PRMC].tail(1))
steady_afterxtal_PRMC_fme_combined = pd.concat(steady_afterxtal_PRMC)
FME_refill = [x+y for x,y in zip(F_xtal_PRMC,F_erupt_PRMC)]
FxtalE_ratio = [fxtale_ratio]*len(F_xtal_PRMC)
FME_df = pd.DataFrame({'xtal=1-F':F_xtal_PRMC,'erupt=E':F_erupt_PRMC,'refill=1-F+E':FME_refill,'(1-F)/E':FxtalE_ratio})  

    
# %% pure fractional crystallization 
magma = LLDinput.loc[index,['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']]
magma['FeO'] = magma['FeO']*0.9  # assume FeO/FeO(t) = 0.9 for MORB
magma_trace = LLDinput.loc[index,['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']]
system_components = magma.to_dict()
T_system_components = oxideToComponent(system_components)
t_start = get_first_T(T_system_components, P, kdCalc = kdCalc_langmuir1992)
t_stop = t_start - 250
tstep = 0.5
trace_components = magma_trace.to_dict()
fl,fa_dict,major_oxide_dict,major_phase_oxide_dict,trace_dict = frac_model_trange(t_start, t_stop, tstep, system_components,trace_components,kd_trace,P,kdCalc = kdCalc_langmuir1992) 

T_df = pd.DataFrame(np.arange(t_start,t_stop,-tstep))
T_df = T_df-273.15
T_df.columns = ['T_C']
fl_dict = {'fl':fl}
fl_df = pd.DataFrame(fl_dict)
fa_df = pd.DataFrame(fa_dict)
major_oxide_df = pd.DataFrame(major_oxide_dict)
major_ol_oxide_df = pd.DataFrame(major_phase_oxide_dict['ol'])
major_cpx_oxide_df = pd.DataFrame(major_phase_oxide_dict['cpx'])
major_plg_oxide_df = pd.DataFrame(major_phase_oxide_dict['plg'])
trace_dict_df = pd.DataFrame(trace_dict)

LLD_df_frac = pd.concat([T_df,fl_df,fa_df,major_oxide_df,major_ol_oxide_df,major_cpx_oxide_df,major_plg_oxide_df,trace_dict_df],axis=1)
LLD_df_frac.columns = ['T_C','f_liq','f_plg','f_cpx','f_ol','liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
               'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','olSiO2',\
                   'olTiO2','olAl2O3','olFeO','olMgO','olK2O','olMnO','olNa2O','olP2O5','olCaO','olNiO',\
                       'cpxSiO2','cpxTiO2','cpxAl2O3','cpxFeO','cpxMgO','cpxK2O','cpxMnO','cpxNa2O',\
                           'cpxP2O5','cpxCaO','cpxNiO','plgSiO2','plgTiO2','plgAl2O3','plgFeO','plgMgO',\
                               'plgK2O','plgMnO','plgNa2O','plgP2O5','plgCaO','plgNiO','liq_La','liq_Ce',\
                                   'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                                       'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
LLD_df_frac['liq_FeOt'] = LLD_df_frac['liq_FeO']/0.9  # assume FeO/FeO(t) = 0.9 for MORB
LLD_df_frac['Fo'] = 100/(1+LLD_df_frac['olFeO']/LLD_df_frac['olMgO']*40.3/71.84)
LLD_df_frac['olNippm'] = LLD_df_frac['olNiO']*58.6934/74.69*10**4
LLD_df_frac['olMnppm'] = LLD_df_frac['olMnO']*54.938/70.94*10**4
LLD_df_frac['liq_Ni'] = LLD_df_frac['liq_NiO']*58.6934/74.69*10**4
LLD_df_frac['liq_FeOtMnO'] = LLD_df_frac['liq_FeOt']/LLD_df_frac['liq_MnO']


# %% pure equilibrium crystallization 
magma = LLDinput.loc[index,['SiO2','TiO2','Al2O3','FeO','MgO','K2O','MnO','Na2O','P2O5','CaO','NiO']]
magma['FeO'] = magma['FeO']*0.9  # assume FeO/FeO(t) = 0.9 for MORB
magma_trace = LLDinput.loc[index,['La','Ce','Nd','Sm','Eu','Gd','Dy','Er','Yb','Lu','Sr','Th','Ba','Rb','Ta','Nb','D=0','Pr','Pb','U']]
system_components = magma.to_dict()
T_system_components = oxideToComponent(system_components)
t_start = get_first_T(T_system_components, P, kdCalc = kdCalc_langmuir1992)
t_stop = t_start -250
tstep = 0.5
trace_components = magma_trace.to_dict()
fl,fa_dict,major_oxide_dict,major_phase_oxide_dict,trace_dict = eq_model_trange(t_start, t_stop,tstep,system_components,trace_components,kd_trace,P,kdCalc = kdCalc_langmuir1992) 

T_df = pd.DataFrame(np.arange(t_start,t_stop,-tstep))
T_df = T_df-273.15
T_df.columns = ['T_C']
fl_dict = {'fl':fl}
fl_df = pd.DataFrame(fl_dict)
fa_df = pd.DataFrame(fa_dict)
major_oxide_df = pd.DataFrame(major_oxide_dict)
major_ol_oxide_df = pd.DataFrame(major_phase_oxide_dict['ol'])
major_cpx_oxide_df = pd.DataFrame(major_phase_oxide_dict['cpx'])
major_plg_oxide_df = pd.DataFrame(major_phase_oxide_dict['plg'])
trace_dict_df = pd.DataFrame(trace_dict)

LLD_df_equ = pd.concat([T_df,fl_df,fa_df,major_oxide_df,major_ol_oxide_df,major_cpx_oxide_df,major_plg_oxide_df,trace_dict_df],axis=1)
LLD_df_equ.columns = ['T_C','f_liq','f_plg','f_cpx','f_ol','liq_SiO2','liq_TiO2','liq_Al2O3','liq_FeO',\
               'liq_MgO','liq_K2O','liq_MnO','liq_Na2O','liq_P2O5','liq_CaO','liq_NiO','olSiO2',\
                   'olTiO2','olAl2O3','olFeO','olMgO','olK2O','olMnO','olNa2O','olP2O5','olCaO','olNiO',\
                       'cpxSiO2','cpxTiO2','cpxAl2O3','cpxFeO','cpxMgO','cpxK2O','cpxMnO','cpxNa2O',\
                           'cpxP2O5','cpxCaO','cpxNiO','plgSiO2','plgTiO2','plgAl2O3','plgFeO','plgMgO',\
                               'plgK2O','plgMnO','plgNa2O','plgP2O5','plgCaO','plgNiO','liq_La','liq_Ce',\
                                   'liq_Nd','liq_Sm','liq_Eu','liq_Gd','liq_Dy','liq_Er','liq_Yb','liq_Lu','liq_Sr',\
                                       'liq_Th','liq_Ba','liq_Rb','liq_Ta','liq_Nb','liq_D=0','liq_Pr','liq_Pb','liq_U']
LLD_df_equ['liq_FeOt'] = LLD_df_equ['liq_FeO']/0.9  # assume FeO/FeO(t) = 0.9 for MORB
LLD_df_equ['Fo'] = 100/(1+LLD_df_equ['olFeO']/LLD_df_equ['olMgO']*40.3/71.84)
LLD_df_equ['olNippm'] = LLD_df_equ['olNiO']*58.6934/74.69*10**4
LLD_df_equ['olMnppm'] = LLD_df_equ['olMnO']*54.938/70.94*10**4
LLD_df_equ['liq_Ni'] = LLD_df_equ['liq_NiO']*58.6934/74.69*10**4
LLD_df_equ['liq_FeOtMnO'] = LLD_df_equ['liq_FeOt']/LLD_df_equ['liq_MnO']

    


