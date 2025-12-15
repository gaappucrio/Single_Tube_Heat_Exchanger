import streamlit as st
from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.animation import FuncAnimation
from seaborn.palettes import blend_palette
import tempfile


def run_simulation(L, r, n, m, Cp, rho, Ti, T0, q_fluxo, t_final, dt):
    dx = L / n
    x = np.linspace(dx/2, L-dx/2, n)
    T = np.ones(n) * T0
    t = np.arange(0, t_final, dt)
    
    # Criando a figura para o gráfico em regime permanente
    fig_permanente = plt.figure(figsize=(8, 6))

    # Função que define a EDO para a variação da temperatura
    def dTdt_function(T, t):
        dTdt = np.zeros(n)
        dTdt[1:n] = (m*Cp*(T[0:n-1]-T[1:n])+q_fluxo*2*np.pi*r*dx)/(rho*Cp*dx*np.pi*r**2)
        dTdt[0] = (m*Cp*(Ti-T[0])+q_fluxo*2*np.pi*r*dx)/(rho*Cp*dx*np.pi*r**2)
        return dTdt
    
    # Resolvendo a EDO usando o odeint
    T_out = odeint(dTdt_function, T, t)
    T_out = T_out
    
    # Criação do DataFrame
    df_Temp = pd.DataFrame(np.array(T_out), columns=x)
    
    # Criando uma paleta de cores
    paleta_calor = blend_palette(['yellow', 'orange','red'], as_cmap=True, n_colors=100)
    
    # Função que atualiza o plot
    def update_plot(t):
        plt.clf()
        line = pd.DataFrame(df_Temp.iloc[t, :]).T
        sns.heatmap(line, cmap=paleta_calor)
        plt.title(f'Time: {t} (s)')
        plt.gca().set_xticklabels(['{:.2f}'.format(val) for val in x])

    #Criando figura para a animação
    fig_animacao = plt.figure(figsize=(8, 6))
    
    # Criando a animação
    ani = FuncAnimation(fig_animacao, update_plot, frames=df_Temp.shape[0], repeat=False)
    # Criar arquivo temporário
    with tempfile.NamedTemporaryFile(suffix=".gif", delete=False) as tmpfile:
        ani.save(tmpfile.name, writer="pillow", fps=10)
        tmpfile.seek(0)
        gif_bytes = tmpfile.read()

    # Salvar a animação como gif
    #save = ani.save('Temperature_Variation_Case_I.gif', writer='pillow', fps=10)
    
    # Exibindo a simulação
    with st.expander("Real-time Simulation Visualization for the Fluid (Click here to view)"):
        st.write('Temperature variation of the fluid passing through the heat exchanger over time and along its length.')
        st.write('Time is shown above the GIF in seconds. Temperatures in Kelvin are represented on the variable scale of the y-axis. The heat exchanger length is shown in meters on the GIF’s x-axis.')
        st.image(gif_bytes, caption="Temperature variation – Case I")

    #Exibindo o gráfico de variação da temperatura ao longo do comprimento em regime permanente
    plt.figure(fig_permanente)
    plt.plot(x, df_Temp.iloc[-1, :], color='blue')  
    plt.xlabel('Length (m)')
    plt.ylabel('Temperature (K)')
    plt.title('Fluid temperature along the exchanger length under steady-state conditions.')
    st.pyplot(plt)

st.title('TROCAL Simulator - Simulation of a single tube heat exchanger')
st.write('This is a simulator of a single tube heat exchanger which heats a fluid as it flows through it. When running the simulation, you will be able to visualize the temperature profile of the fluid você poderá visualizar o perfil de temperatura do fluido along the heat exchanger as time progresses. You can also view the steady-state temperatures along the length of the heat exchanger.')
st.write('Below is an illustrative figure of this heat exchanger, created by the authors.')
st.image('Case 1.png', use_column_width=True)
st.write('One application of this configuration is the use of coils in on-demand (instantaneous) heating systems. These coils consist of tubes or tube assemblies through which fluids flow and are heated by an external heat source.')
st.write('This case can also represent any industrial tubulation applied for heating fluids through an external heating source.')
st.write('This simulator employs the following energy balance equation for the fluid flowing through the heat exchanger, based on the principle of energy conservation:')
st.image('Equation Case 1.jpg', use_column_width=True)
st.write('ATENTION: At the bottom of this page, you will also find a button that runs the simulation using a predefined example (‘Run standard example’). This example takes approximately 30 seconds to complete, depending on your connection speed. If you prefer to use your own input values, select the ‘Run simulation’ button. It is recommended to use between 10 and 30 nodes, depending on the specific case being analyzed.”')


st.title('Simulation Input Parameters')
# Valores input
L = st.number_input('Tube length (m)', min_value=0.0)
r = st.number_input('Tube radius (m)', min_value=0.0)
n = st.number_input('Number of nodes for discretization', min_value=1)
m = st.number_input('Mass Flow (kg/s)', min_value=0.0)
Cp = st.number_input('Specific heat capacity of the fluid (J/kg.K)', min_value=0.0)
rho = st.number_input('Specific mass of the fluid (kg/m³)', min_value=0.0)
Ti = st.number_input('Inlet fluid temperature (K)')
T0 = st.number_input('Initial temperature of the heat exchanger (K)')
q_fluxo = st.number_input('Heat flux (W/m²)', min_value=0.0)
t_final = st.number_input('Simulation time (s)', min_value=0.0)
dt = st.number_input('Time step (s)', min_value=0.0)

if st.button('Run Simulation'):
    run_simulation(L, r, n, m, Cp, rho, Ti, T0, q_fluxo, t_final, dt)
elif st.button('Run standard example'):
    run_simulation(10, 0.1, 10, 3, 4180, 995.61, 400, 300, 10000, 210, 1)




