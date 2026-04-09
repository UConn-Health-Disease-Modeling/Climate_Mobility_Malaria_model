#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# ode_model.py
import numpy as np

def malaria_model(df_sim, climate_series, mobility, params):
    
    Temperature_data_urban, Temperature_data_rural = climate_series["T_u"], climate_series["T_r"];
    Rainfall_data_urban, Rainfall_data_rural = climate_series["R_u"], climate_series["R_r"];
    U_to_R1_series, U_to_R2_value = mobility["m_UR1_t"], mobility["m_UR2"];
    R_to_U1_series, R_to_U2_value = mobility["m_RU1_t"], mobility["m_RU2"];
    
    for i in range(1, len(df_sim['tempo'])):
        Tt_U = Temperature_data_urban[i-1]
        R_U = Rainfall_data_urban[i-1]
        Tt_R = Temperature_data_rural[i-1]
        R_R = Rainfall_data_rural[i-1]

        params['R_to_U1'] = R_to_U1_series[i]
        params['R_to_U2'] = R_to_U2_value

        params['U_to_R1'] = U_to_R1_series[i]
        params['U_to_R2'] = U_to_R2_value
        
        if i > 22*7:
            R_U = Rainfall_data_urban[i-22*7]
            R_R = Rainfall_data_rural[i-22*7]
        
        if Temperature_data_urban[i] < 1.0:
            Tt_U = Temperature_data_urban[i-1]
        if Temperature_data_rural[i] < 1.0:
            Tt_R = Temperature_data_rural[i-1]
            
        if i > 2*7:
            Tt_U = Temperature_data_urban[i-2*7]
            Tt_R = Temperature_data_rural[i-2*7]

        # Urban Mosquito birth rate as a function of temperature
        pE_U = (4*params['pME']/params['RLE_U']**2)*R_U*(params['RLE_U'] - R_U)
        pLR_U = (4*params['pML']/params['RLE_U']**2)*R_U*(params['RLE_U'] - R_U)
        pP_U = (4*params['pMP']/params['RLE_U']**2)*R_U*(params['RLE_U'] - R_U)

        # Rural Mosquito birth rate as a function of temperature
        pE_R = (4*params['pME']/params['RLE_R']**2)*R_R*(params['RLE_R'] - R_R)
        pLR_R = (4*params['pML']/params['RLE_R']**2)*R_R*(params['RLE_R'] - R_R)
        pP_R = (4*params['pMP']/params['RLE_R']**2)*R_R*(params['RLE_R'] - R_R)
        
        if R_U > params['RLE_U']:
            pE_U = 0.0
            pLR_U = 0.0
            pP_U = 0.0
        if R_R > params['RLE_R']:
            pE_R = 0.0
            pLR_R = 0.0
            pP_R = 0.0
        
        pLT_U = np.exp(-(params['alpham']*Tt_U + params['betam']))
        tauL_U = 1/(params['alpham']*Tt_U + params['betam'])
        a_U = params['c_u1']*(0.000203*(Tt_U**2 - 11.7*Tt_U)*(np.sqrt(42.3 - Tt_U)))
        um_U = params['c_u3']*(-np.log(-0.000828*Tt_U**2 + 0.0367*Tt_U + 0.522))
        #Pc = 1e6
        BE_U = (-0.153*Tt_U**2 + 8.61*Tt_U - 97.7)/um_U
        tau_U = 1/(-0.00094*Tt_U**2 + 0.049*Tt_U - 0.0552)
        deltam_U = params['c_u5']*(BE_U*pE_U*pLR_U*pLT_U*pP_U/tau_U)
        
        pLT_R = np.exp(-(params['alpham']*Tt_R + params['betam']))
        tauL_R = 1/(params['alpham']*Tt_R + params['betam'])
        a_R = params['c_r1']*(0.000203*(Tt_R**2 - 11.7*Tt_R)*(np.sqrt(42.3 - Tt_R)))
        um_R = params['c_r2']*(-np.log(-0.000828*Tt_R**2 + 0.0367*Tt_R + 0.522))
        #Pc = 1e6
        BE_R = (-0.153*Tt_R**2 + 8.61*Tt_R - 97.7)/um_R
        tau_R = 1/(-0.00094*Tt_R**2 + 0.049*Tt_R - 0.0552)
        deltam_R = params['c_r3']*(BE_R*pE_R*pLR_R*pLT_R*pP_R/tau_R)


        # Variables
        Sh_UU = df_sim['Sh_UU'][i-1]
        Sh_UR = df_sim['Sh_UR'][i-1]
        Ih_UU = df_sim['Ih_UU'][i-1]
        Ih_UR = df_sim['Ih_UR'][i-1]
        L_UU = df_sim['L_UU'][i-1]
        L_UR = df_sim['L_UR'][i-1]
        P_UU = df_sim['P_UU'][i-1]
        P_UR = df_sim['P_UR'][i-1]
        Sm_U = df_sim['Sm_U'][i-1]
        Im_U = df_sim['Im_U'][i-1]
        #Aw_U = df_sim['Aw_U'][i-1]
        Nh_U = df_sim['Nh_U'][i-1]
        Nm_U = df_sim['Nm_U'][i-1]
        
        Sh_RR = df_sim['Sh_RR'][i-1]
        Sh_RU = df_sim['Sh_RU'][i-1]
        Ih_RR = df_sim['Ih_RR'][i-1]
        Ih_RU = df_sim['Ih_RU'][i-1]
        L_RR = df_sim['L_RR'][i-1]
        L_RU = df_sim['L_RU'][i-1]
        P_RR = df_sim['P_RR'][i-1]
        P_RU = df_sim['P_RU'][i-1]
        Sm_R = df_sim['Sm_R'][i-1]
        Im_R = df_sim['Im_R'][i-1]
        #Aw_R = df_sim['Aw_R'][i-1]
        Nh_R = df_sim['Nh_R'][i-1]
        Nm_R = df_sim['Nm_R'][i-1]

        # Human presence in each region
        Nh_U = Sh_UU + Ih_UU + L_UU + P_UU + Sh_RU + Ih_RU + L_RU + P_RU
        Nh_R = Sh_RR + Ih_RR + L_RR + P_RR + Sh_UR + Ih_UR + L_UR + P_UR

        Nm_U = Sm_U + Im_U
        m_U = Nm_U/Nh_U

        Nm_R = Sm_R + Im_R
        m_R = Nm_R/Nh_R

        # Adjust force of infection using residence time
        lambda_Uh = m_U * a_U * params['b_U'] * Im_U / Nm_U
        lambda_Rh = m_R * a_R * params['b_R'] * Im_R / Nm_R
        
        # Infected human presence
        Ih_U_pres = Ih_UU + Ih_RU
        Ih_R_pres = Ih_RR + Ih_UR
        P_U_pres = P_UU + P_RU
        P_R_pres = P_RR + P_UR

        # Infection by Ih to mosquitoes
        lambda_Uv_Ih = (a_U * params['cs'] * params['sigmav'] + a_U * params['ca'] * (1 - params['sigmav'])) * Ih_U_pres / Nh_U
        lambda_Rv_Ih = (a_R * params['cs'] * params['sigmav'] + a_R * params['ca'] * (1 - params['sigmav'])) * Ih_R_pres / Nh_R
        lambda_Uv_P = a_U * params['cs'] * params['epsv'] / params['kappa'] * (1 - params['phi']) * P_U_pres / Nh_U
        lambda_Rv_P = a_R * params['cs'] * params['epsv'] / params['kappa'] * (1 - params['phi']) * P_R_pres / Nh_R

        # Calculating variations for Urban
        dSh_UU = -lambda_Uh*Sh_UU + (1 - params['c_u4']*params['eta_U']*params['sigmav'])*(1 - params['phiu'])*params['rv']*Ih_UU + params['uvL']*L_UU + (1 - params['phit']*(1 - params['phi']))/params['kappa']*P_UU - params['U_to_R1']*Sh_UU + params['U_to_R2']*Sh_UR
        dSh_UR = -params['c_u2_1']*lambda_Rh*Sh_UR + (1 - params['c_u4']*params['eta_U']*params['sigmav'])*(1 - params['phiu'])*params['rv']*Ih_UR + params['uvL']*L_UR + (1 - params['phit']*(1 - params['phi']))/params['kappa']*P_UR + params['U_to_R1']*Sh_UU - params['U_to_R2']*Sh_UR
        dIh_UU = lambda_Uh*Sh_UU - (1 - params['c_u4']*params['eta_U']*params['sigmav'])*params['rv']*Ih_UU - params['c_u4']*params['eta_U']*params['sigmav']*params['gammav']*Ih_UU + params['psi']*L_UU + lambda_Uh*L_UU - params['U_to_R1']*Ih_UU + params['U_to_R2']*Ih_UR
        dIh_UR = params['c_u2_1']*lambda_Rh*Sh_UR - (1 - params['c_u4']*params['eta_U']*params['sigmav'])*params['rv']*Ih_UR - params['c_u4']*params['eta_U']*params['sigmav']*params['gammav']*Ih_UR + params['psi']*L_UR + params['c_u2_1']*lambda_Rh*L_UR + params['U_to_R1']*Ih_UU - params['U_to_R2']*Ih_UR
        dCases_U = lambda_Uh*Sh_UU + params['psi']*L_UU + lambda_Uh*L_UU + lambda_Rh*Sh_UR + params['psi']*L_UR + lambda_Rh*L_UR
        dL_UU = (1 - params['c_u4']*params['eta_U']*params['sigmav'])*params['phiu']*params['rv']*Ih_UU + params['phit']*(1 - params['phi'])/params['kappa']*P_UU - params['uvL']*L_UU - params['psi']*L_UU - lambda_Uh*L_UU - params['U_to_R1']*L_UU + params['U_to_R2']*L_UR
        dL_UR = (1 - params['c_u4']*params['eta_U']*params['sigmav'])*params['phiu']*params['rv']*Ih_UR + params['phit']*(1 - params['phi'])/params['kappa']*P_UR - params['uvL']*L_UR - params['psi']*L_UR - params['c_u2_1']*lambda_Rh*L_UR + params['U_to_R1']*L_UU - params['U_to_R2']*L_UR
        dP_UU = params['c_u4']*params['eta_U']*params['sigmav']*params['gammav']*Ih_UU - P_UU/params['kappa'] - params['U_to_R1']*P_UU + params['U_to_R2']*P_UR
        dP_UR = params['c_u4']*params['eta_U']*params['sigmav']*params['gammav']*Ih_UR - P_UR/params['kappa'] + params['U_to_R1']*P_UU - params['U_to_R2']*P_UR
        dSm_U = -um_U*Sm_U + deltam_U - lambda_Uv_Ih*Sm_U - lambda_Uv_P*Sm_U
        dIm_U = -um_U*Im_U + lambda_Uv_Ih*Sm_U + lambda_Uv_P*Sm_U
        #dAw_U = params['omega']*(params['ini_awa_U'] + params['xi']*(Ih_UU + Ih_UR))*(1 - Aw_U) - params['lambda']*Aw_U

        # Calculating variations for Rural
        dSh_RR = -lambda_Rh*Sh_RR + (1 - params['c_r4']*params['eta_R']*params['sigmav'])*(1 - params['phiu'])*params['rv']*Ih_RR + params['uvL']*L_RR + (1 - params['phit']*(1 - params['phi']))/params['kappa']*P_RR - params['R_to_U1']*Sh_RR + params['R_to_U2']*Sh_RU
        dSh_RU = -params['c_u2']*lambda_Uh*Sh_RU + (1 - params['c_r5']*params['eta_R']*params['sigmav'])*(1 - params['phiu'])*params['rv']*Ih_RU + params['uvL']*L_RU + (1 - params['phit']*(1 - params['phi']))/params['kappa']*P_RU + params['R_to_U1']*Sh_RR - params['R_to_U2']*Sh_RU
        dIh_RR = lambda_Rh*Sh_RR - (1 - params['c_r4']*params['eta_R']*params['sigmav'])*params['rv']*Ih_RR - params['c_r4']*params['eta_R']*params['sigmav']*params['gammav']*Ih_RR + params['psi']*L_RR + lambda_Rh*L_RR - params['R_to_U1']*Ih_RR + params['R_to_U2']*Ih_RU
        dIh_RU = params['c_u2']*lambda_Uh*Sh_RU - (1 - params['c_r5']*params['eta_R']*params['sigmav'])*params['rv']*Ih_RU - params['c_r5']*params['eta_R']*params['sigmav']*params['gammav']*Ih_RU + params['psi']*L_RU + params['c_u2']*lambda_Uh*L_RU + params['R_to_U1']*Ih_RR - params['R_to_U2']*Ih_RU
        dCases_R = lambda_Rh*Sh_RR + params['psi']*L_RR + lambda_Rh*L_RR + params['c_u2']*lambda_Uh*Sh_RU + params['psi']*L_RU + params['c_u2']*lambda_Uh*L_RU
        dL_RR = (1 - params['c_r4']*params['eta_R']*params['sigmav'])*params['phiu']*params['rv']*Ih_RR + params['phit']*(1 - params['phi'])/params['kappa']*P_RR - params['uvL']*L_RR - params['psi']*L_RR - lambda_Rh*L_RR - params['R_to_U1']*L_RR + params['R_to_U2']*L_RU
        dL_RU = (1 - params['c_r5']*params['eta_R']*params['sigmav'])*params['phiu']*params['rv']*Ih_RU + params['phit']*(1 - params['phi'])/params['kappa']*P_RU - params['uvL']*L_RU - params['psi']*L_RU - params['c_u2']*lambda_Uh*L_RU + params['R_to_U1']*L_RR - params['R_to_U2']*L_RU
        dP_RR = params['c_r4']*params['eta_R']*params['sigmav']*params['gammav']*Ih_RR - P_RR/params['kappa'] - params['R_to_U1']*P_RR + params['R_to_U2']*P_RU
        dP_RU = params['c_r5']*params['eta_R']*params['sigmav']*params['gammav']*Ih_RU - P_RU/params['kappa'] + params['R_to_U1']*P_RR - params['R_to_U2']*P_RU
        dSm_R = -um_R*Sm_R + deltam_R - lambda_Rv_Ih*Sm_R - lambda_Rv_P*Sm_R
        dIm_R = -um_R*Im_R + lambda_Rv_Ih*Sm_R + lambda_Rv_P*Sm_R
        #dAw_R = params['omega']*(params['ini_awa_R'] + params['xi']*(Ih_RR + Ih_RU))*(1 - Aw_R) - params['lambda']*Aw_R

        # Applying variations to Urban
        df_sim.at[i, 'Sh_UU'] = Sh_UU + dSh_UU
        df_sim.at[i, 'Sh_UR'] = Sh_UR + dSh_UR
        df_sim.at[i, 'Ih_UU'] = Ih_UU + dIh_UU
        df_sim.at[i, 'Ih_UR'] = Ih_UR + dIh_UR
        df_sim.at[i, 'L_UU'] = L_UU + dL_UU
        df_sim.at[i, 'L_UR'] = L_UR + dL_UR
        df_sim.at[i, 'P_UU'] = P_UU + dP_UU
        df_sim.at[i, 'P_UR'] = P_UR + dP_UR
        #df_sim.at[i, 'Aw_U'] = Aw_U + dAw_U
        df_sim.at[i, 'Sm_U'] = Sm_U + dSm_U
        df_sim.at[i, 'Im_U'] = Im_U + dIm_U
        df_sim.at[i, 'Cases_U'] = dCases_U
        df_sim.at[i, 'Nm_U'] = Nm_U
        df_sim.at[i, 'Nh_U'] = Nh_U
        m_U = df_sim.at[i, 'Nm_U']/df_sim.at[i, 'Nh_U']

        # Applying variations to Rural
        df_sim.at[i, 'Sh_RR'] = Sh_RR + dSh_RR
        df_sim.at[i, 'Sh_RU'] = Sh_RU + dSh_RU
        df_sim.at[i, 'Ih_RR'] = Ih_RR + dIh_RR
        df_sim.at[i, 'Ih_RU'] = Ih_RU + dIh_RU
        df_sim.at[i, 'L_RR'] = L_RR + dL_RR
        df_sim.at[i, 'L_RU'] = L_RU + dL_RU
        df_sim.at[i, 'P_RR'] = P_RR + dP_RR
        df_sim.at[i, 'P_RU'] = P_RU + dP_RU
        #df_sim.at[i, 'Aw_R'] = Aw_R + dAw_R
        df_sim.at[i, 'Sm_R'] = Sm_R + dSm_R
        df_sim.at[i, 'Im_R'] = Im_R + dIm_R
        df_sim.at[i, 'Cases_R'] = dCases_R
        df_sim.at[i, 'Nm_R'] = Nm_R
        df_sim.at[i, 'Nh_R'] = Nh_R
        m_R = df_sim.at[i, 'Nm_R']/df_sim.at[i, 'Nh_R']

    return df_sim

