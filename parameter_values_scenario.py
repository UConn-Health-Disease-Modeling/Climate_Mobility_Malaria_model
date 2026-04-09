#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Dictionary of parameters used in malaria_model, grouped by category
params = {
    # Climate control parameters
        'pME': 0.90,          
        'pML': 0.25,          
        'pMP': 0.75,          
        'RLE_U': 53.19243,
        'RLE_R': 33.0,    
        'alpham': 0.0554,     
        'betam': -0.06737,     
    # Resistance model parameters
        'm_U': 2435/625,
        'm_R': 2435/625,
        'b_U': 0.005095394,   
        'b_R': 0.005595394,   
        'eta_U': 0.01,        
        'eta_R': 0.05,        
        'sigmav': 0.33,       
        'rv': 1.0/60.0,       
        'gammav': 1.0/9.0,    
        'phi': 0.95,
        'kappa': 30.0,        
        'epsv': 2.1,          
        'ca': 0.12,           
        'cs': 0.4,            
        'psi': 1.0/60.0,      
        'phit': 0.21,         
        'phiu': 0.4,          
        'uvL': 1.0/425.0,    
    
    # Intervention control params
        "c_u1": 1,
        "c_u2": 1, 
        "c_u2_1": 1, 
        "c_u3": 1, 
        "c_u4": 1, 
        "c_u5": 1, 
        "c_r1": 1, 
        "c_r2": 1, 
        "c_r3": 1, 
        "c_r4":1, 
        "c_r5":1
}







