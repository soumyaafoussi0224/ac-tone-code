

import math
import csv
import numpy as np
import matplotlib.pyplot as plt




FLUID = "Acetone"

# Feed conditions
T_FEED_C = 10.0
P_FEED_ATM = 1.0
M_FEED_KG_H = 30.0

# Splitter
SPLIT_1 = 0.40
SPLIT_2 = 0.60

# Heating pressure
P_HEATER_ATM = 1.0

# Parametric study
N_POINTS = 25
PRESSURE_STEP_ATM = 0.1
P_START_ATM = 0.1

# Heater efficiency
HEATER_EFFICIENCY = 1.0

# Acetone molecular weight
MW_ACETONE = 58.080  # kg/kmol




ANTOINE_A = 4.41877
ANTOINE_B = 1312.253
ANTOINE_C = -32.445


# Approximate liquid heat capacity of acetone
# around the temperature range of interest.
# Unit: kJ/(kg.K)

CP_LIQUID = 2.15


# Approximate latent heat of vaporization at normal boiling point
# Unit: kJ/kg

HVAP_NORMAL = 518.0


# Normal boiling point of acetone
T_BOILING_NORMAL_K = 329.15


# ============================================================
# 3. UNIT CONVERSION FUNCTIONS
# ============================================================

def c_to_k(T_C):
    return T_C + 273.15


def k_to_c(T_K):
    return T_K - 273.15


def atm_to_pa(P_atm):
    return P_atm * 101325.0


def pa_to_atm(P_pa):
    return P_pa / 101325.0


def kg_h_to_kg_s(m):
    return m / 3600.0


def kg_h_to_kmol_h(m):
    return m / MW_ACETONE


def kJ_h_to_kW(Q):
    return Q / 3600.0


# ============================================================
# 4. ANTOINE EQUATION
# ============================================================

def vapor_pressure_atm(T_K):
    """
    Saturation pressure of acetone in atm.
    """

    log10_P = ANTOINE_A - (
        ANTOINE_B / (T_K + ANTOINE_C)
    )

    return 10.0 ** log10_P


def boiling_temperature_K(P_atm):
    """
    Calculate acetone boiling temperature
    at a given pressure using Antoine equation.
    """

    # Antoine:
    #
    # log10(P) = A - B/(T+C)
    #
    # Rearrangement:
    #
    # T = B/(A-log10(P)) - C

    T_K = (
        ANTOINE_B /
        (ANTOINE_A - math.log10(P_atm))
        - ANTOINE_C
    )

    return T_K


def boiling_temperature_C(P_atm):
    return k_to_c(
        boiling_temperature_K(P_atm)
    )


# ============================================================
# 5. LATENT HEAT MODEL
# ============================================================

def latent_heat_kJ_kg(T_K):
    """
    Approximate latent heat using a Watson-type correlation.

    Hvap(T) = Hvap(Tb) *
              [(1-Tr)/(1-Trb)]^0.38

    Tc acetone ≈ 508.1 K
    """

    T_CRITICAL_K = 508.10

    Tr = T_K / T_CRITICAL_K
    Trb = T_BOILING_NORMAL_K / T_CRITICAL_K

    if Tr >= 1.0:
        return 0.0

    exponent = 0.38

    Hvap = HVAP_NORMAL * (
        (1.0 - Tr) /
        (1.0 - Trb)
    ) ** exponent

    return Hvap


# ============================================================
# 6. SENSIBLE HEATING
# ============================================================

def sensible_heat_kJ_kg(T_initial_C, T_final_C):
    """
    Energy required to heat liquid acetone
    from T_initial to T_final.

    Q = Cp * DeltaT
    """

    delta_T = T_final_C - T_initial_C

    return CP_LIQUID * delta_T


# ============================================================
# 7. COMPLETE VAPORIZATION ENERGY
# ============================================================

def vaporization_energy(
        m_kg_h,
        T_initial_C,
        P_atm):

    # Boiling temperature
    T_boil_C = boiling_temperature_C(P_atm)

    T_boil_K = c_to_k(T_boil_C)

    # Sensible heating
    q_sensible = sensible_heat_kJ_kg(
        T_initial_C,
        T_boil_C
    )

    # Latent heat
    q_latent = latent_heat_kJ_kg(
        T_boil_K
    )

    # Total specific energy
    q_total = q_sensible + q_latent

    # Total energy per hour
    Q_kJ_h = m_kg_h * q_total

    # Thermal power
    Q_kW = kJ_h_to_kW(Q_kJ_h)

    # Correct for heater efficiency
    Q_required_kW = Q_kW / HEATER_EFFICIENCY

    return {
        "T_boiling_C": T_boil_C,
        "T_boiling_K": T_boil_K,
        "q_sensible_kJ_kg": q_sensible,
        "q_latent_kJ_kg": q_latent,
        "q_total_kJ_kg": q_total,
        "Q_kJ_h": Q_kJ_h,
        "Q_kW": Q_kW,
        "Q_required_kW": Q_required_kW
    }


# ============================================================
# 8. SPLITTER CALCULATION
# ============================================================

def splitter(m_feed):

    m1 = m_feed * SPLIT_1
    m2 = m_feed * SPLIT_2

    return m1, m2


# ============================================================
# 9. STREAM PROPERTY CALCULATION
# ============================================================

def stream_information(
        name,
        T_C,
        P_atm,
        mass_flow):

    T_K = c_to_k(T_C)

    molar_flow = kg_h_to_kmol_h(mass_flow)

    Psat = vapor_pressure_atm(T_K)

    boiling_T = boiling_temperature_C(P_atm)

    # Simple phase identification
    if T_C < boiling_T - 0.01:
        phase = "LIQUID"
    elif abs(T_C - boiling_T) <= 0.01:
        phase = "SATURATED / EQUILIBRIUM"
    else:
        phase = "VAPOR"

    return {
        "Name": name,
        "Temperature_C": T_C,
        "Temperature_K": T_K,
        "Pressure_atm": P_atm,
        "Pressure_kPa": atm_to_pa(P_atm) / 1000,
        "MassFlow_kg_h": mass_flow,
        "MassFlow_kg_s": kg_h_to_kg_s(mass_flow),
        "MolarFlow_kmol_h": molar_flow,
        "MolarFlow_mol_s": molar_flow * 1000 / 3600,
        "MolecularWeight_kg_kmol": MW_ACETONE,
        "VaporPressure_atm": Psat,
        "BoilingTemperature_C": boiling_T,
        "Phase": phase
    }





def run_process():

    print("\n")
    print("=" * 70)
    print("          ACETONE PROCESS SIMULATION")
    print("=" * 70)


    feed = stream_information(
        "FEED",
        T_FEED_C,
        P_FEED_ATM,
        M_FEED_KG_H
    )

  

    m1, m2 = splitter(M_FEED_KG_H)

    stream_1 = stream_information(
        "STREAM-1 / 40%",
        T_FEED_C,
        P_FEED_ATM,
        m1
    )

    stream_2_in = stream_information(
        "STREAM-2 / 60% INLET",
        T_FEED_C,
        P_HEATER_ATM,
        m2
    )

    # --------------------------------------------------------
    # HEATER
    # --------------------------------------------------------

    energy = vaporization_energy(
        m2,
        T_FEED_C,
        P_HEATER_ATM
    )

    stream_2_out = stream_information(
        "VAPOR OUTLET",
        energy["T_boiling_C"],
        P_HEATER_ATM,
        m2
    )

    # --------------------------------------------------------
    # DISPLAY FEED
    # --------------------------------------------------------

    print("\nFEED STREAM")
    print("-" * 70)

    for key, value in feed.items():
        print(f"{key:35s}: {value}")

    # --------------------------------------------------------
    # SPLITTER
    # --------------------------------------------------------

    print("\nSPLITTER")
    print("-" * 70)

    print(
        f"Feed mass flow        = {M_FEED_KG_H:.4f} kg/h"
    )

    print(
        f"Outlet 1 (40%)        = {m1:.4f} kg/h"
    )

    print(
        f"Outlet 2 (60%)        = {m2:.4f} kg/h"
    )

    print(
        f"Mass balance error     = "
        f"{M_FEED_KG_H - m1 - m2:.6e} kg/h"
    )

    
    print("\nBOILING POINT AT 1 ATM")
    print("-" * 70)

    print(
        f"Pressure               = {P_HEATER_ATM:.4f} atm"
    )

    print(
        f"Boiling temperature    = "
        f"{energy['T_boiling_C']:.4f} °C"
    )

    print(
        f"Boiling temperature    = "
        f"{energy['T_boiling_K']:.4f} K"
    )

   
    print("\nHEATER")
    print("-" * 70)

    print(
        f"Mass flow              = {m2:.4f} kg/h"
    )

    print(
        f"Initial temperature    = {T_FEED_C:.4f} °C"
    )

    print(
        f"Final temperature      = "
        f"{energy['T_boiling_C']:.4f} °C"
    )

    print(
        f"Sensible heat          = "
        f"{energy['q_sensible_kJ_kg']:.4f} kJ/kg"
    )

    print(
        f"Latent heat            = "
        f"{energy['q_latent_kJ_kg']:.4f} kJ/kg"
    )

    print(
        f"Total specific heat    = "
        f"{energy['q_total_kJ_kg']:.4f} kJ/kg"
    )

    print(
        f"Energy required        = "
        f"{energy['Q_kJ_h']:.4f} kJ/h"
    )

    print(
        f"Thermal power          = "
        f"{energy['Q_kW']:.4f} kW"
    )

    print(
        f"Required heater power  = "
        f"{energy['Q_required_kW']:.4f} kW"
    )


    print("\nVAPOR OUTLET")
    print("-" * 70)

    for key, value in stream_2_out.items():
        print(f"{key:35s}: {value}")

    return {
        "feed": feed,
        "stream1": stream_1,
        "stream2_in": stream_2_in,
        "stream2_out": stream_2_out,
        "energy": energy
    }




def pressure_sensitivity():

    """
    25 pressure points.

    Starting pressure = 0.1 atm
    Increment          = 0.1 atm

    Therefore:
    0.1, 0.2, 0.3, ..., 2.5 atm
    """

    pressures = []

    boiling_temperatures = []
    sensible_energies = []
    latent_energies = []
    total_energies = []

    m_vapor = M_FEED_KG_H * SPLIT_2

    for i in range(N_POINTS):

        P = P_START_ATM + i * PRESSURE_STEP_ATM

        T_boil_C = boiling_temperature_C(P)

        T_boil_K = c_to_k(T_boil_C)

        q_sensible = sensible_heat_kJ_kg(
            T_FEED_C,
            T_boil_C
        )

        q_latent = latent_heat_kJ_kg(
            T_boil_K
        )

        q_total = q_sensible + q_latent

        pressures.append(P)
        boiling_temperatures.append(T_boil_C)
        sensible_energies.append(q_sensible)
        latent_energies.append(q_latent)
        total_energies.append(q_total)

    return (
        np.array(pressures),
        np.array(boiling_temperatures),
        np.array(sensible_energies),
        np.array(latent_energies),
        np.array(total_energies)
    )




def display_pressure_table():

    (
        pressures,
        temperatures,
        sensible,
        latent,
        total
    ) = pressure_sensitivity()

    print("\n")
    print("=" * 90)
    print("PRESSURE / BOILING TEMPERATURE STUDY")
    print("=" * 90)

    print(
        f"{'Point':>5s}"
        f"{'P (atm)':>12s}"
        f"{'P (kPa)':>12s}"
        f"{'Tb (°C)':>14s}"
        f"{'Qsens':>14s}"
        f"{'Qlatent':>14s}"
        f"{'Qtotal':>14s}"
    )

    print("-" * 90)

    for i in range(len(pressures)):

        print(
            f"{i+1:5d}"
            f"{pressures[i]:12.3f}"
            f"{pressures[i] * 101.325:12.3f}"
            f"{temperatures[i]:14.3f}"
            f"{sensible[i]:14.3f}"
            f"{latent[i]:14.3f}"
            f"{total[i]:14.3f}"
        )




def plot_results():

    (
        pressures,
        temperatures,
        sensible,
        latent,
        total
    ) = pressure_sensitivity()

    

    plt.figure(figsize=(9, 6))

    plt.plot(
        temperatures,
        pressures,
        marker="o"
    )

    plt.xlabel(
        "Boiling Temperature (°C)"
    )

    plt.ylabel(
        "Pressure (atm)"
    )

    plt.title(
        "Acetone Vapor Pressure / Boiling Curve"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.show()

   

    plt.figure(figsize=(9, 6))

    plt.plot(
        pressures,
        temperatures,
        marker="o"
    )

    plt.xlabel(
        "Pressure (atm)"
    )

    plt.ylabel(
        "Boiling Temperature (°C)"
    )

    plt.title(
        "Acetone Boiling Temperature vs Pressure"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.show()

    

    plt.figure(figsize=(9, 6))

    plt.plot(
        pressures,
        sensible,
        marker="o",
        label="Sensible heat"
    )

    plt.plot(
        pressures,
        latent,
        marker="s",
        label="Latent heat"
    )

    plt.plot(
        pressures,
        total,
        marker="^",
        label="Total heat"
    )

    plt.xlabel(
        "Pressure (atm)"
    )

    plt.ylabel(
        "Specific Energy (kJ/kg)"
    )

    plt.title(
        "Energy Requirement vs Pressure"
    )

    plt.legend()

    plt.grid(True)

    plt.tight_layout()

    plt.show()




def export_csv(filename="acetone_pressure_study.csv"):

    (
        pressures,
        temperatures,
        sensible,
        latent,
        total
    ) = pressure_sensitivity()

    with open(
        filename,
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.writer(file)

        writer.writerow([
            "Point",
            "Pressure_atm",
            "Pressure_kPa",
            "Boiling_Temperature_C",
            "Sensible_Heat_kJ_kg",
            "Latent_Heat_kJ_kg",
            "Total_Heat_kJ_kg"
        ])

        for i in range(len(pressures)):

            writer.writerow([
                i + 1,
                pressures[i],
                pressures[i] * 101.325,
                temperatures[i],
                sensible[i],
                latent[i],
                total[i]
            ])

    print(
        f"\nResults exported to: {filename}"
    )



def display_parameters():

    print("\n")
    print("=" * 70)
    print("SIMULATION PARAMETERS")
    print("=" * 70)

    parameters = {
        "Fluid": FLUID,
        "Feed temperature (°C)": T_FEED_C,
        "Feed pressure (atm)": P_FEED_ATM,
        "Feed mass flow (kg/h)": M_FEED_KG_H,
        "Splitter outlet 1 (%)": SPLIT_1 * 100,
        "Splitter outlet 2 (%)": SPLIT_2 * 100,
        "Heater pressure (atm)": P_HEATER_ATM,
        "Heater efficiency (%)": HEATER_EFFICIENCY * 100,
        "Number of pressure points": N_POINTS,
        "Pressure starting point (atm)": P_START_ATM,
        "Pressure increment (atm)": PRESSURE_STEP_ATM,
        "Molecular weight (kg/kmol)": MW_ACETONE,
        "Liquid Cp (kJ/kg.K)": CP_LIQUID,
        "Normal latent heat (kJ/kg)": HVAP_NORMAL
    }

    for name, value in parameters.items():

        print(
            f"{name:40s}: {value}"
        )



if __name__ == "__main__":

    display_parameters()

    results = run_process()

    display_pressure_table()

    plot_results()

    export_csv()

    print("\n")
    print("=" * 70)
    print("SIMULATION FINISHED")
    print("=" * 70)