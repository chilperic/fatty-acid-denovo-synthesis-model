#!/usr/bin/env python
# coding: utf-8
"""
This model differs from the other as the CoA inhibit the binding of acetyl-CoA and malonyl-CoA 
in the  competitve manner.
"""

import numpy as np
import pandas as pd
import scipy.integrate as sp
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import matplotlib.patches as mpatch
import matplotlib.lines as lines
import scipy.stats
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.lines as lines
import scipy.stats
from math import erf, log10, ceil
from collections.abc import Callable
import sys
from typing import Tuple
from math import erf, log10
from scipy.interpolate import lagrange
from scipy.linalg import null_space

####for fitting ##################

import random as rd

plt.style.use("default")
np.warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

get_ipython().run_line_magic("matplotlib", "inline")


# In[2]:


compounds = [
    "AcetylCoA",
    "MalonylCoA",
    # "C4FA",
    # "C6FA",
    # "C8FA",
    # "C10FA",
    # "C12FA",
    "C14FA",
    "C16FA",
    "C18FA",
    "CoA",
    "E",
    "EC2",
    "EC4",
    "EC6",
    "EC8",
    "EC10",
    "EC12",
    "EC14",
    "EC16",
    "EC18",
    "ECoA",
    "NADPH",
]

n_compounds = len(compounds)
n_compounds


# In[3]:


color = [
    ["b", " "],
    ["c", "-"],
    ["xkcd:cranberry", "-"],
    ["m", "-"],
    ["xkcd:sepia", "-"],
    ["xkcd:very dark brown", "-"],
    ["xkcd:tomato", "-"],
    ["tab:olive", "--"],
    ["xkcd:azure", "-"],
    ["xkcd:chartreuse", "-"],
    ["xkcd:cyan", "--"],
    ["xkcd:fuchsia", "--"],
    ["xkcd:beige", "--"],
    ["xkcd:gold", "-."],
    ["xkcd:plum", "--"],
    ["xkcd:red", "-"],
    ["xkcd:coral", "--"],
    ["xkcd:orange", "--"],
    ["xkcd:black", "--"],
    ["xkcd:tan", "--"],
    ["xkcd:khaki", "--"],
    ["xkcd:green", "-"],
    ["xkcd:sienna", "-."],
    ["xkcd:navy", "--"],
    ["xkcd:salmon", "--"],
    ["xkcd:goldenrod", "--"],
]


# color = ["b", "k", "g", "r", "c", "m", "y"]


# ## Generic functions to be used within the model

# In[4]:


def mass_actions(k, *args):
    """
    this function simulate mass action kinetics
    """
    expr = k
    for arg in args:
        expr = expr * arg
    return expr


def reversible_mass_action(kf, kb, list_subtrates, list_products):
    forward = kf
    back_ward = kb
    for sub in list_subtrates:
        forward = forward * sub
    for prod in list_products:
        back_ward = back_ward * prod
    expr = forward - back_ward
    return expr
    

def rev_mass_action(kf, kb, s1,s2,p1,p2=1):
    expr = kf*s1*s2 - kb*p1*p2
    return expr


def mass_actions_with_inhibition(k, ki, S):
    expr = k * S / (1 + ki * S)
    return expr


def mass_action_inhib(k, ki, s1, s2, s3,  i):
    "In the case there is only one substrate, set s3 =1"
    expr = k*s1*s2*s3/(1 +ki*i)
    return expr

def constant_influx(k):
    """
    This function simulate the zero oder
    kinetics constant influx
    """
    return k


def proportional_outflux(k, S):
    """This function is just like mass action, 
    it simulate proportional outfluxes"""
    return k * S


def ord_bibi_uniuni_pingpong(
    Vf, Vr, Kma, Kmb, Kmc, Kmp, Kmq, Kmr, Keq,
):
    pass


def logistic_rate(r, K, x):
    rate = r * x * (1 - x / K)
    return rate


def michaelis(Vmax, Km, S):
    rate = (Vmax * S) / (S + Km)
    return rate


def normal(n, mean=16, std=0.5):
    expr = scipy.stats.norm.pdf(n, mean, std)
    return expr


def skewed_normal(
    n, c, mean=16, std=5,
):
    expr = 2 * normal(n, mean, std) * normal(c * n, mean, std)
    return expr


def derived_rate_release(n, mean=16, std=0.5, t=1):
    r = -(1 / t) * np.log(1 - normal(n, mean, std))
    return r


def derived_rate_attaching(n, mean=4, std=5, t=1):
    r = -(1 / t) * np.log(1 - normal(n, mean, std))
    return r


def exp_rate_release(n, k, m):
    expr = m * np.exp(k * n)
    return expr


def exp_rate_attaching(n, k, m):
    expr = m * np.exp(-k * n)
    return expr


def solubility(n, a=370.89 * 1e-6, b=1.385):
    expr = a * np.exp(-b * n)
    return expr


def conservation(total, *args):
    expr = total - sum(args)
    return expr


conservation(1, 0.1, 0.1, 0.1, 0.1)


# ## Utilitity functions

# In[5]:


def simulate_time_course(model, par, y0, max_time, ncp=0, ncp_list=None):
    """
    time is a numpy.ndarray that can be generated with numpy.linspace(lower, upper, number_of_points)
    y0 = initial conditions (numpy.ndarray)  with the dimension of the number of metabolites in the compounds 
    par =  is a dictionnary of parameters
    color: a list of list defining the color for the plot  (see color above)

    """
    problem = Explicit_Problem(lambda t, y: model(t, y, par), y0)
    integrator = CVode(problem)
    integrator.rtol = 1e-12
    integrator.atol = 1e-12 
    integrator.verbosity = 50  # 10 (loud) - 50 (quiet)
    if ncp_list is None:
        t, y = integrator.simulate(max_time, ncp)
    else:
        max_time = ncp_list[-1]
        t, y = integrator.simulate(max_time, ncp, ncp_list)

    return t, y


def plot_time_course(time, y, compounds, color):
    """ y here is the solution of the ODE system
    """
    n_compounds = len(compounds)

    # plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(10, 7))
    for i in range(n_compounds):
        plt.plot(time, y[:, i], color=color[i][0], linewidth=3)
    plt.xlabel("Time (a.u.) ", fontsize=12)
    plt.ylabel("Concemtration (a.u.) ", fontsize=12)
    plt.title("Concentrations over time ", fontsize=14)
    ax.legend(compounds, loc="upper left", bbox_to_anchor=(1.1, 1), borderaxespad=0)
    plt.tight_layout()
    plt.savefig("Concentrations.png", bbox_inches="tight")
    plt.show()


def check_steady_state(y, tol=1e-7):
    """
    function to check if the system have reached the
    steady state by comparing the norm of differences in the last values
     y: numpy.array
     tol: float
    """

    if np.linalg.norm(y[-1] - y[-2], ord=2) < tol:
        return True
    return False


def simulate_steady_state_scipy(model, y0, par, dt=1e-6, error=1e-8, max_steps=1e6):
    r = sp.ode(model)
    r.set_integrator("lsoda", method="bdf")
    r.set_f_params(par)
    r.set_initial_value(y0, 0.0)

    y1 = y0.copy()
    nsteps = 0
    while r.successful and nsteps < max_steps:
        y2 = r.integrate(r.t + dt)
        if np.linalg.norm(y1 - y2, ord=2) < error:
            return y2
        else:
            y1 = y2
        nsteps += 1
    raise ValueError("Maximal simulation steps exceeded")


def simulate_steady_state_assimulo(
    model, y0, par, n_iteration=3, time_init=100, tol=1e-7
):
    """
    This function return last value if the system reaches the steady state, otherwise it return None
    Requires Explicit_Problem and CVode from assimulo package

    """
    problem = Explicit_Problem(lambda t, y: model(t, y, par), y0)
    integrator = CVode(problem)
    integrator.verbosity = 50  # 10 (loud) - 50 (quiet)
    integrator.rtol = 10e-14
    integrator.atol = 10e-14 * np.ones(n_compounds)

    time_to_simulate = time_init
    for i in range(n_iteration):
        try:
            t, y = integrator.simulate(time_to_simulate)
            if check_steady_state(y, tol):
                return y[-1]
            time_to_simulate = 10 * time_to_simulate
        except CVodeError:
            return None
    return None


def parameter_scan(rhs, y0, par, par_name, par_values, s2ss_function):
    """
    This function change the value of the parameters and record the
    steady state values and the rates values at the steady state.
    It returns two dictionnaries namely the dictionnary of steady state concentrations
    where the dict_keys are the values of the varied parameter and
    the dictionnary of steady state rates
    where the dict_keys are the values of the varied parameter
    and the dict_values are steady state rates
    rhs = is the model (system of ODE)
    y0  =  initial condition for the model
    par = dictionnary of parameters
    par_name =  parameter to vary
    par_values =  list or 1D numpy.array of the values in which the parameter will be varried
    """
    par = par.copy()
    steady_state_values = {}
    steady_state_fluxes = {}
    for value in par_values:
        par[par_name] = value
        y_steady_state = s2ss_function(rhs, y0, par)
        if y_steady_state is not None:
            steady_state_values[value] = y_steady_state
            steady_state_fluxes[value] = rates(y_steady_state, par)
    return (
        pd.DataFrame(steady_state_values, index=compounds).T,
        pd.DataFrame(steady_state_fluxes).T,
    )


def get_subplots(
    number_of_plots: int,
    number_of_columns: int = 1,
    base_size_row: int = 3,
    base_size_column: int = 4,
) -> Tuple[plt.Figure, np.ndarray]:

    """
    This function helps the define properly the axis for the plots
    number_of_plots: int, (That is the number of plots we want)
    number_of_columns: int (= 1,)
    fig_size_row: int = 12,
    fig_size_column: int = 8,
    """
    nrows = ceil(number_of_plots / number_of_columns)
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=number_of_columns,
        figsize=(base_size_column * number_of_columns, base_size_row * nrows),
        squeeze=False,
    )
    axs = axs.flatten()
    print(axs)
    return fig, axs


def plot_groups(
    data: pd.DataFrame,
    groups: list,
    ref_value: float = None,
    ncols: int = 2,
    base_size_row: int = 3,
    base_size_column: int = 4,
):
    fig, axs = get_subplots(
        len(groups),
        number_of_columns=ncols,
        base_size_column=base_size_column,
        base_size_row=base_size_row,
    )
    for ax, group in zip(axs, groups):
        ax.plot(data.loc[:, group], "-o", linewidth=3)
        if ref_value is not None:
            ax.axvline(x=ref_value, c="r")
            group.append("Ref_value")

        ax.legend(
            group,
            bbox_to_anchor=(1.01, 1),
            loc="upper left",
            borderaxespad=0,
            fontsize=14,
        )
    fig.tight_layout()
    return fig, axs


# In[ ]:


# ## <div align="center"> The different rate equations
#
# In this set up we supposed that all the reactions the simple mass action kinetics with. We also assume that the total amount of $ATP$ is fixed as well as the total amount
#
# 1. We still have to handle the $CO2/ HCO3$
# 2.  Same with ACC1 and FAS
#
#
# One way for the enzymes will be to include the release of the free enzyme at the end of the reactions.

# In[6]:


def rates(y, par):
    (
        AcetylCoA,
        MalonylCoA,
        # C4FACoA,
        # C6FACoA,
        # C8FACoA,
        # C10FACoA,
        # C12FACoA,
        C14FACoA,
        C16FACoA,
        C18FACoA,
        CoA,
        E,
        EC2,
        EC4,
        EC6,
        EC8,
        EC10,
        EC12,
        EC14,
        EC16,
        EC18,
        ECoA,
        NADPH,
    ) = y
    
    # kd1 = par["nu"] * derived_rate_release(4, std =par["std"])
    # kd2 = par["nu"] * derived_rate_release(6, std =par["std"])
    # kd3 = par["nu"] * derived_rate_release(8, std =par["std"])
    # kd4 = par["nu"] * derived_rate_release(10, std =par["std"])
    # kd5 = par["nu"] * derived_rate_release(12, std =par["std"])
    # kd6 = par["nu"] * derived_rate_release(14, std =par["std"])
    # kd7 = par["nu"] * derived_rate_release(16, std =par["std"])
    # kd8 = par["nu"] * derived_rate_release(18, std =par["std"])

    rate = {
        "AcetylCoA__EC2":mass_actions(par["beta"], E, AcetylCoA),
        "EC2__EC4": mass_actions(par["kappa"], EC2, NADPH,MalonylCoA),
        "EC4__EC6":  mass_actions(par["kappa"] ,EC4, NADPH, MalonylCoA),
        "EC6__EC8":  mass_actions(par["kappa"] ,EC6, NADPH, MalonylCoA),
        "EC8__EC10":  mass_actions(par["kappa"] ,EC8, NADPH, MalonylCoA),
        "EC10__EC12": mass_actions(par["kappa"] ,EC10,NADPH, MalonylCoA),
        "EC12__EC14":  mass_actions(par["kappa"] ,EC12,NADPH, MalonylCoA),
        "EC14__EC16":  mass_actions(par["kappa"] ,EC14,NADPH, MalonylCoA),
        "EC16__EC18":  mass_actions(par["delta"] ,EC16,NADPH, MalonylCoA),
        ##########################################################################
        # "EC4__C4FA": mass_actions(par["kd1"], EC4),
        # "EC6__C6FA": mass_actions(par["kd2"], EC6),
        # "EC8__C8FA": mass_actions(par["kd3"], EC8),
        # "EC10__C10FA": mass_actions(par["kd4"], EC10),
        # "EC12__C12FA": mass_actions(par["kd5"], EC12),
        "EC14__C14FA": mass_actions(par["kd6"], EC14),
        "EC16__C16FA": mass_actions(par["kd7"], EC16),
        "EC18__C18FA": mass_actions(par["kd8"], EC18),
        ######################################################################
        # "C4_FA__EC4": mass_actions(par["ka1"], E, C4FACoA),
        # "C6_FA__EC6": mass_actions(par["ka2"], E, C6FACoA),
        # "C8_FA__EC8": mass_actions(par["ka3"], E, C8FACoA),
        # "C10_FA__EC10": mass_actions(par["ka4"], E, C10FACoA),
        # "C12_FA__EC12": mass_actions(par["ka5"], E, C12FACoA),
        # "C14_FA__EC14": mass_actions(par["ka6"], E, C14FACoA),
        # "C16_FA__EC16": mass_actions(par["ka7"], E, C16FACoA),
        # "C18_FA__EC18": mass_actions(par["ka8"], E, C18FACoA),
        #########################################################################
        "E__ECoA": rev_mass_action(kf =par["kf_CoA"], kb =par["kb_CoA"], s1 =E,s2 =CoA,p1 =ECoA,p2=1),
    }
    # print(E)

    return rate


# # <div align="center"> The ODE system

# In[ ]:


# In[ ]:


# In[7]:


def system(time, y, par):
    get_rates = rates(par=par, y=y)
    derivatives = {
        "dAcetylCoA": -get_rates["AcetylCoA__EC2"],
        "dMalonylCoA": -(
            get_rates["EC2__EC4"]
            + get_rates["EC4__EC6"]
            + get_rates["EC6__EC8"]
            + get_rates["EC8__EC10"]
            + get_rates["EC10__EC12"]
            + get_rates["EC12__EC14"]
            + get_rates["EC14__EC16"]
            + get_rates["EC16__EC18"]
        ),
        # "dC4FA": get_rates["EC4__C4FA"] ,
        # "dC6FA": get_rates["EC6__C6FA"] ,
        # "dC8FA": get_rates["EC8__C8FA"] ,
        # "dC10FA": get_rates["EC10__C10FA"],
        # "dC12FA": get_rates["EC12__C12FA"],
        "dC14FA": get_rates["EC14__C14FA"],
        "dC16FA": get_rates["EC16__C16FA"],
        "dC18FA": get_rates["EC18__C18FA"],

        "dCoA": -get_rates["E__ECoA"]
        +(get_rates["AcetylCoA__EC2"]
          + get_rates["EC2__EC4"] 
          + get_rates["EC4__EC6"]
          + get_rates["EC6__EC8"]
          + get_rates["EC8__EC10"]
          + get_rates["EC10__EC12"]
          + get_rates["EC12__EC14"]
          + get_rates["EC14__EC16"]
          + get_rates["EC16__EC18"]),
        "dE": -get_rates["AcetylCoA__EC2"] - get_rates["E__ECoA"]+
     (
            # get_rates["EC4__C4FA"]
            # + get_rates["EC6__C6FA"]
            # + get_rates["EC8__C8FA"]
            # + get_rates["EC10__C10FA"]
            # + get_rates["EC12__C12FA"]
            + get_rates["EC14__C14FA"]
            + get_rates["EC16__C16FA"]
            + get_rates["EC18__C18FA"]
        ),
        
        "dEC2": get_rates["AcetylCoA__EC2"] - get_rates["EC2__EC4"],
        "dEC4": get_rates["EC2__EC4"]
        - get_rates["EC4__EC6"],
        "dEC6": get_rates["EC4__EC6"]
        - get_rates["EC6__EC8"],
        "dEC8": get_rates["EC6__EC8"]
        - get_rates["EC8__EC10"],
        "dEC10": get_rates["EC8__EC10"]
        - get_rates["EC10__EC12"],
        "dEC12": get_rates["EC10__EC12"]
        - get_rates["EC12__EC14"],
        "dEC14": get_rates["EC12__EC14"]
        - get_rates["EC14__EC16"]
        - get_rates["EC14__C14FA"],
        "dEC16": get_rates["EC14__EC16"]
        - get_rates["EC16__EC18"]
        - get_rates["EC16__C16FA"],
        "dEC18": get_rates["EC16__EC18"]
        - get_rates["EC18__C18FA"],
        "dECoA": get_rates["E__ECoA"],
        "dNADPH": -2
        * (
            get_rates["EC2__EC4"]
            + get_rates["EC4__EC6"]
            + get_rates["EC6__EC8"]
            + get_rates["EC8__EC10"]
            + get_rates["EC10__EC12"]
            + get_rates["EC12__EC14"]
            + get_rates["EC14__EC16"]
            + get_rates["EC16__EC18"]
        ),
    }
    # da = [derivatives[f"dEC{i}"] for i in range(2, 20, 2)]
    # da.append(derivatives["dECoA"])
    # print(sum(da)+ derivatives["dE"])

    return np.fromiter(derivatives.values(), dtype=float)

# In[8]:


def fatty_acid_of_interest(system, par, y0, max_time, ncp_list=None):
    time, y = simulate_time_course(
        model=system, par=par, y0=y0, max_time=max_time, ncp_list=ncp_list
    )
    y = y.T
    expr = np.array([y[2], y[3], y[4]]).T
    return expr


def fatty_acid_of_interest2(system, par, y0, time):
    y = sp.odeint(func=system, y0=y0, t=time, args=(par,), tfirst=True)
    y = y.T
    expr = np.array([y[2], y[3], y[4]]).T
    return expr


def total_fatty_acid(system, par, y0, max_time, ncp_list=None):
    time, y = simulate_time_course(
        model=system, par=par, y0=y0, max_time=max_time, ncp_list=ncp_list
    )
    y = y.T
    expr = np.array(sum((y[2], y[3], y[4]))).T
    return expr


def goodness_of_fit(data, model):
    """This function measure the goodness of the 
    fit between the data and the simulated model
    """
    expr = np.sqrt(np.mean(np.square(np.subtract(data, model))))
    return expr
    
    
    
    
def specific_fatty_acid(
    specific_compound: str,
    compounds: list[str],
    system,
    par: dict[str, float],
    y0: float,
    max_time: float,
    ncp_list=None,
):
    time, y = simulate_time_course(
        model=system, par=par, y0=y0, max_time=max_time, ncp_list=ncp_list
    )
    y = y.T
    y_dict = dict(zip(compounds, y))
    expr = y_dict[specific_compound]

    return expr

    expr = np.array(sum((y[2], y[3], y[4]))).T
    return expr




def time_course_production(
    var_substrate: str,
    prod: str,
    compounds: list[str],
    system,
    y0_dict: dict[str, float],
    par: dict[str, float],
    lower_val: float,
    upper_val: float,
    number_of_values: int,
    max_time: float,
    ncp_list: list,
):
    """
    The function   help to build a family of curve of time couse 
    by varying a specific substrates

    Returns:
        _type_: _description_
    """
    y0_dict1 = y0_dict.copy()
    list_substrate_values = np.linspace(lower_val, upper_val, number_of_values)
    data_dict = {}
    fig, ax = plt.subplots(figsize=(8, 6))
    for i in list_substrate_values:
        y0_dict1[var_substrate] = i

        z0 = np.array(list(y0_dict1.values()))

        time_course_prod = specific_fatty_acid(
            prod, compounds, system, par, z0, max_time, ncp_list
        )
        data_dict[f"{i}"] = time_course_prod

        ax.plot(
            ncp_list, data_dict[f"{i}"], "-o", label=var_substrate + " =" + f"{i}",
        )
    
    ax.set_xlabel("time (s)")
    ax.set_ylabel(f"concentration of {prod} ($\mu M$)")
    ax.set_title(f" Time course of {prod} at different concentration of  {var_substrate} ")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 0.5))

    return data_dict
    
    
    
def steady_state_kinetics(
    var_substrate: str,
    prod: str,
    compounds: list[str],
    system,
    y0_dict: dict[str, float],
    par: dict[str, float],
    lower_val: float,
    upper_val: float,
    number_of_values: int,
    max_time: float,
    ncp_list: list,
):
    """
    The function   help to build  setady state kinetics of a particular product 
    base onn the variation of a particular substrate

    Returns:
        _type_: _description_
    """
    y0_dict1 = y0_dict.copy()
    list_substrate_values = np.linspace(lower_val, upper_val, number_of_values)
    # data_dict = {}
    rate_kinetics = []
    Delta_t = ncp_list[-1] - ncp_list[0]
    for i in list_substrate_values:
        y0_dict1[var_substrate] = i

        z0 = np.array(list(y0_dict1.values()))

        time_course_prod = specific_fatty_acid(
            prod, compounds, system, par, z0, max_time, ncp_list
        )
        # data_dict[f"{i}"] = time_course_prod

        Delta_p = time_course_prod[-1] - time_course_prod[0]
        rate_of_change = Delta_p / Delta_t
        rate_kinetics.append(rate_of_change)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(list_substrate_values, rate_kinetics)
    ax.plot(list_substrate_values, rate_kinetics)
    ax.set_xlabel(f"{var_substrate} ($\mu M$)")
    ax.set_ylabel(f"rate kinetics of {prod} ($\mu M \cdot s^{-1}$)")
    ax.set_title(f"Steady state kinetics of {prod} as function of {var_substrate} ")

    return ax
    
    



# %%
