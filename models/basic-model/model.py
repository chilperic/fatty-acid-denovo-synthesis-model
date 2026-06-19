#!/usr/bin/env python
# coding: utf-8



"""
This model does not include the CoA inhibition or any other types of inhibition .
"""

from __future__ import annotations

from math import ceil
from typing import Callable, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate as sp
import scipy.stats
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode
from assimulo.solvers.sundials import CVodeError


get_ipython().run_line_magic('matplotlib', 'inline')




compounds = [
    "AcetylCoA",
    "MalonylCoA",
    "C4FA",
    "C6FA",
    "C8FA",
    "C10FA",
    "C12FA",
    "C14FA",
    "C16FA",
    "C18FA",
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
    "NADPH",
]

n_compounds = len(compounds)
n_compounds



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
    ["xkcd:sienna", '-.'],
    ["xkcd:navy", "--"],
    ["xkcd:salmon", "--"],
    ["xkcd:goldenrod", "--"],
]



def mass_actions(k, *args):
    """
    this function simulate mass action kinetics
    """
    expr = k 
    for arg in args:
        expr = expr *arg
    return expr


def reversible_mass_action(kf, kb,list_subtrates, list_products):
    forward = kf
    back_ward = kb
    for sub in list_subtrates:
        forward = forward*sub
    for prod in list_products:
        back_ward = back_ward*prod
    expr =  forward -back_ward
    return expr
    
    
    
def mass_actions_with_inhibition(k,ki, S):
    expr = k*S/(1+ki*S)
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
    Vf,
    Vr,
    Kma,
    Kmb,
    Kmc,
    Kmp,
    Kmq,
    Kmr,
    Keq,
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
    n,
    c,
    mean=16,
    std=5,
):
    expr = 2 * normal(n, mean, std) * normal(c * n, mean, std)
    return expr


def derived_rate_release(n, mean=16, std=0.5, t=1):
    r = -(1 / t) * np.log(1 - normal(n, mean, std))
    return r


def derived_rate_attaching(n, mean=4, std=5, t=1):
    r = -(1 / t) * np.log(1 - normal(n, mean, std))
    return r


def exp_rate_release(n, k,m):
    expr = m*np.exp(k*n)
    return expr


def exp_rate_attaching(n, k,m):
    expr = m*np.exp(-k * n)
    return expr

def solubility(n, a =370.89*1e-6, b= 1.385):
    expr = a*np.exp(-b*n)
    return expr

def conservation(total, *args):
    expr =  total -sum(args)
    return expr


conservation(1, 0.1, 0.1, 0.1, 0.1)





def simulate_time_course(model, par, y0, max_time, ncp=0, ncp_list =None):
    """
    time is a numpy.ndarray that can be generated with numpy.linspace(lower, upper, number_of_points)
    y0 = initial conditions (numpy.ndarray)  with the dimension of the number of metabolites in the compounds 
    par =  is a dictionnary of parameters
    color: a list of list defining the color for the plot  (see color above)

    """
    problem = Explicit_Problem(lambda t, y: model(t, y, par), y0)
    integrator = CVode(problem)
    integrator.rtol =1e-15
    integrator.atol = 1e-15
    integrator.verbosity = 50  # 10 (loud) - 50 (quiet)
    if ncp_list is None:
        t, y = integrator.simulate(max_time, ncp)
    else:
        max_time = ncp_list[-1]
        t, y = integrator.simulate(max_time,ncp, ncp_list)
   
    return t,y



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
    ax.legend(compounds,
              loc="upper left",
              bbox_to_anchor=(1.1, 1),
              borderaxespad=0)
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
    integrator.rtol =10e-10
    integrator.atol = 10e-10

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
    ref_value:float =None,
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
            
        ax.legend(group, bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0, fontsize=14)
    fig.tight_layout()
    return fig, axs    



def rates(y, par):
    (
        AcetylCoA,
        MalonylCoA,
        C4FACyt,
        C6FACyt,
        C8FACyt,
        C10FACyt,
        C12FACyt,
        C14FACyt,
        C16FACyt,
        C18FACyt,
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
        NADPH,
    ) = y
   

    rate = {
        
        "AcetylCoA__EC2": mass_actions(par["beta"], E, AcetylCoA),
        "EC2__EC4": mass_actions(par["kappa"], EC2, NADPH, MalonylCoA),
        "EC4__EC6": mass_actions(par["kappa"], EC4, NADPH, MalonylCoA),
        "EC6__EC8": mass_actions(par["kappa"], EC6, NADPH, MalonylCoA),
        "EC8__EC10": mass_actions(par["kappa"], EC8, NADPH, MalonylCoA),
        "EC10__EC12": mass_actions(par["kappa"], EC10, NADPH, MalonylCoA),
        "EC12__EC14": mass_actions(par["kappa"], EC12, NADPH, MalonylCoA),
        "EC14__EC16": mass_actions(par["kappa"], EC14, NADPH, MalonylCoA),
        "EC16__EC18": mass_actions(par["delta"], EC16, NADPH, MalonylCoA),
        ##########################################################################
        "EC4__C4FA": mass_actions(par["kd1"], EC4),
        "EC6__C6FA": mass_actions(par["kd2"], EC6),
        "EC8__C8FA": mass_actions(par["kd3"], EC8),
        "EC10__C10FA": mass_actions(par["kd4"], EC10),
        "EC12__C12FA": mass_actions(par["kd5"], EC12),
        "EC14__C14FA": mass_actions(par["kd6"], EC14),
        "EC16__C16FA": mass_actions(par["kd7"], EC16),
        "EC18__C18FA": mass_actions(par["kd8"], EC18),
        ######################################################################
    
    }
 

    return rate





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
        "dC4FA": get_rates["EC4__C4FA"] ,
        "dC6FA": get_rates["EC6__C6FA"] ,
        "dC8FA": get_rates["EC8__C8FA"] ,
        "dC10FA": get_rates["EC10__C10FA"] ,
        "dC12FA": get_rates["EC12__C12FA"] ,
        "dC14FA": get_rates["EC14__C14FA"] ,
        "dC16FA": get_rates["EC16__C16FA"] ,
        "dC18FA": get_rates["EC18__C18FA"] ,
        "dE": -get_rates["AcetylCoA__EC2"] + (get_rates["EC4__C4FA"]
                +get_rates["EC6__C6FA"] 
                +get_rates["EC8__C8FA"] 
                +get_rates["EC10__C10FA"]
                +get_rates["EC12__C12FA"]
                +get_rates["EC14__C14FA"]
                +get_rates["EC16__C16FA"]
                +get_rates["EC18__C18FA"]),
        "dEC2": get_rates["AcetylCoA__EC2"] - get_rates["EC2__EC4"],
        "dEC4": get_rates["EC2__EC4"]
        - get_rates["EC4__EC6"]
        - get_rates["EC4__C4FA"],
        "dEC6": get_rates["EC4__EC6"]
        - get_rates["EC6__EC8"]
        - get_rates["EC6__C6FA"],
        "dEC8": get_rates["EC6__EC8"]
        - get_rates["EC8__EC10"]
        - get_rates["EC8__C8FA"],
        "dEC10": get_rates["EC8__EC10"]
        - get_rates["EC10__EC12"]
        - get_rates["EC10__C10FA"],
        "dEC12": get_rates["EC10__EC12"]
        - get_rates["EC12__EC14"]
        - get_rates["EC12__C12FA"],
        "dEC14": get_rates["EC12__EC14"]
        - get_rates["EC14__EC16"]
        - get_rates["EC14__C14FA"],
        "dEC16": get_rates["EC14__EC16"]
        - get_rates["EC16__EC18"]
        - get_rates["EC16__C16FA"],
        "dEC18": get_rates["EC16__EC18"]
        - get_rates["EC18__C18FA"],
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
  

    return np.fromiter(derivatives.values(), dtype=float)




def fatty_acid_of_interest(system, par, y0, max_time, ncp_list = None):
    time, y = simulate_time_course(model=system, par=par, y0=y0, max_time=max_time, ncp_list = ncp_list)
    y = y.T
    expr = np.array([y[7], y[8], y[9]]).T
    return expr

def fatty_acid_of_interest2(system, par, y0, time):
    y = sp.odeint(func=system, y0=y0, t=time, args=(par,), tfirst=True)
    y = y.T
    expr = np.array([y[7], y[8], y[9]]).T
    return expr


def total_fatty_acid(system, par, y0, max_time, ncp_list = None):
    time, y = simulate_time_course(model=system, par=par, y0=y0, max_time=max_time, ncp_list = ncp_list)
    y = y.T
    expr = np.array(sum((y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9]))).T
    return expr


def goodness_of_fit(data, model):
    """This function measure the goodness of the 
    fit between the data and the simulated model
    """
    expr = np.sqrt(np.mean(np.square(np.subtract(data, model))))
    return expr
