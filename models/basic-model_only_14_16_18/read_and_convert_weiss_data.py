from __future__ import annotations
from pathlib import Path

import pandas as pd


r"""Data from Weiss et al. 2003
 This data present time course data for 240 min  phase for C16:0 production using a purified fatty acid synthase enzyme. 
  
 In this section we get data from this study and convert them to be used for fitting and validate the present
 model. 
 
 To convert everything in the $\mu M$, given that the data are reported in $DPM$, we use the following 
 1. The data on the NEFA are in $DPM$
 
 The convertion is done using the following formula 
 \begin{align}
  data(\mu M) &= \dfrac{data(DPM)}{681\pm 17}
 \end{align}
"""



def read_and_convert_data(
    initial_phase_data_path: Path = Path("../../data/denovo_data_weiss.csv")
    
):

    time_course_data = pd.read_csv("../../data/denovo_data_weiss.csv")
    time_course_data_lwb = time_course_data.copy()
    time_course_data_upb = time_course_data.copy()
    time_course_data_mean = time_course_data.copy()
    time_course_data


    #Convert the data in micro mol. The data is reported in pmol/micro gram FASN. 
    # We would like to have it in $\mu M$ ($\mu mol L^{-1}$)

    def convert_data_in_micro_M(data,conversion):
        concentration_in_micro_Moles =  data/conversion
        return concentration_in_micro_Moles


    conversion_factor_lwb = 681+17

    for i in time_course_data_lwb.columns:
        if i != "time(s)":
            time_course_data_lwb[i] = convert_data_in_micro_M(time_course_data[i],conversion_factor_lwb)
    time_course_data_lwb


    conversion_factor_upb = 681-17

    for i in time_course_data_upb.columns:
        if i != "time(s)":
            time_course_data_upb[i] = convert_data_in_micro_M(time_course_data[i],conversion_factor_upb)
    time_course_data_upb

    
    conversion_factor_mean = 681

    for i in time_course_data_mean.columns:
        if i != "time(s)":
            time_course_data_mean[i] = convert_data_in_micro_M(time_course_data[i],conversion_factor_mean)
    time_course_data_mean



    ############### Correct for initial conditions for  the case where the initial conditions match with the data from the paper 
    deviation_C16 = 1.30
    deviation_C14 = deviation_C16 *1.5/100
    deviation_C18 = deviation_C16 *7/100
    
    time_course_data_mean["C14:0"][0] =deviation_C14
    time_course_data_mean["C16:0"][0] =deviation_C16
    time_course_data_mean["C18:0"][0] =deviation_C18
    time_course_data_mean["C14:0"] =time_course_data_mean["C14:0"]-deviation_C14
    time_course_data_mean["C16:0"] =time_course_data_mean["C16:0"]-deviation_C16
    time_course_data_mean["C18:0"] =time_course_data_mean["C18:0"]-deviation_C18


    time_course_data_mean
    
    time_course_data_lwb_numpy = time_course_data_lwb.to_numpy()
    time_course_data_upb_numpy = time_course_data_upb.to_numpy()
    time_course_data_mean_numpy = time_course_data_mean.to_numpy()
    time_experiment = time_course_data_lwb_numpy[:,0]
    return (
        time_experiment,  # array
        time_course_data_mean,  # Dataframe
       time_course_data_mean_numpy,  # array
    )

   



