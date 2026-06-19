from __future__ import annotations

r"""Data from Topolska et al. 2020

In this section we get data from this study and convert them to be used for fitting and validate the present
model.
To convert everything in the $\mu M$, given that the data are reported in $p mol/ \mu  g FASN$, we use the following
1. The data on the NEFA are in $pmol/ \mu g FASN$
The convertion is done using the following formula
\begin{eqnarray}
    data(\mu M) &=&  \dfrac{data(pmol/ \mu g\cdot FASN)\times 10^{-6}\times massFASN (\mu g)}{vol (L)}\\
    &=& \dfrac{data(pmol/ \mu g\cdot FASN)\times 10^{-6}\times 16.5 \mu g FASN}{200\times 10^{-6}}\\
    &=& 0.0825\times data(pmol/ \mu g\cdot FASN)
\end{eqnarray}
To have data in $\mu M$ , we do the followimg conversion
$$data (\mu M) =   0.0825\times data(pmol/ \mu g\cdot FASN)$$
"""


from pathlib import Path

import pandas as pd


def convert_data_in_micro_M(data, conversion):
    concentration_in_micro_Moles = data * conversion
    return concentration_in_micro_Moles


def read_and_convert_data(
    initial_phase_data_path: Path = Path("../../data/merge_data.csv"),
    tota_fatty_acids_path: Path = Path("../../data/total_fatty acids.csv"),
):
    # Data of C14, C16 and C18 fatty acids produced at the initial phase (120 secs)
    initial_phase_data = pd.read_csv(initial_phase_data_path)

    # Convert the data in micro mol. The data is reported in pmol/micro gram FASN.
    # We would like to have it in $\mu M$ ($\mu mol L^{-1}$)
    conversion_factor = 0.0825

    for i in initial_phase_data.columns:
        if i != "time(s)":
            initial_phase_data[i] = convert_data_in_micro_M(
                initial_phase_data[i], conversion_factor
            )

    initial_phase_data = initial_phase_data.to_numpy()

    total_nefas = pd.read_csv(tota_fatty_acids_path)
    for i in total_nefas.columns:
        if i != "time(s)":
            total_nefas[i] = convert_data_in_micro_M(total_nefas[i], conversion_factor)

    # Time points at which the measurement was done
    time_nefas1_simulation = list(total_nefas["time(s)"].to_numpy())
    total_nefas1 = total_nefas.to_numpy()

    total_nefas1 = total_nefas1[:, 1:]
    total_nefas1 = total_nefas1.flatten()
    return (
        initial_phase_data,  # array
        time_nefas1_simulation,  # array
        total_nefas1,  # array
    )
