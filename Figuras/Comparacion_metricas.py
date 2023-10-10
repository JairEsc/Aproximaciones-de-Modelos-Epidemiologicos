import time
start_time = time.time()
np.random.seed(0)
N=10**3
n=int(N*0.1)
I_0=int(0.01*N)
ksum=1
while(ksum%2!=0):#Que defina una grafica.
        k_dist=np.random.poisson(n,size=N)
        ksum=sum(k_dist)
G=nx.configuration_model(k_dist)
G = nx.Graph(G)#Remover multi-aristas
G.remove_edges_from(nx.selfloop_edges(G))#Remover bucles
pos=nx.spring_layout(G,seed=2,k=3)
print("--- %s Segundos para crear la gráfica ---" % (time.time() - start_time))
X_I_indexes=np.random.choice(N,size=I_0,replace=False).tolist()
for k in range(N):
    if k in X_I_indexes:
        G.nodes[k]['Estado:']='Infeccioso'
    else:
        G.nodes[k]['Estado:']='Susceptible'
print("Empieza")
a_I=1.5
a_R=1.2
import math
beta=75
gamma=1
lammbda_I=((beta/N)*math.gamma(1+1/a_I))
lammbda_R=(gamma*math.gamma(1+1/a_R))
def psi_survival_infection(x):
    return np.exp(-(lammbda_I*x)**(a_I))
def psi_survival_recovery(x):
    return np.exp(-(lammbda_R*x)**(a_R))
def psi_infection(x):
    return (lammbda_I*a_I)*(lammbda_I*x)**(a_I-1)*np.exp(-(lammbda_I*x)**(a_I))
def psi_recovery(x):
    return (lammbda_R*a_R)*(lammbda_R*x)**(a_R-1)*np.exp(-(lammbda_R*x)**(a_R))
times_general=[]
times_general_approx=[]
times_weibull=[]
times_markovian=[]
matriz_trayectorias=[['MG']]
picos={"gen":[],"genApprox":[],"weib":[],"mark":[]}
for k in range(5):
    #------------------------------------ General
    iter_time = time.time()
    H=[]
    while(len(H)<100):
      G_copy=G.copy()
      X_I_indexes_copy=X_I_indexes.copy()
      H=Gillespie_Direct_Method_Non_markovian_network(N=N,X_I_indexes=X_I_indexes_copy,method='general',params=[],densities=[psi_infection,psi_recovery],survivals=[psi_survival_infection,psi_survival_recovery],G_network_initial=G_copy,T=14)
    if(k==0):
        plt.plot(H['time'],H['I(t)'], 'green',linewidth=0.7,label='General (MG)')
    else:
        plt.plot(H['time'],H['I(t)'], 'green',linewidth=0.7)
    i=np.argmax(H['I(t)'])
    picos['gen'].append([H['time'][i],H['I(t)'][i]])
    times_general.append((time.time() - iter_time))
    matriz_trayectorias.append(H['time'])
    matriz_trayectorias.append(H['I(t)'])
    matriz_trayectorias.append(['MG-T'])
    #------------------------------------
    #------------------------------------ General aproximado
    iter_time = time.time()
    H=[]
    while(len(H)<100):
      G_copy=G.copy()
      X_I_indexes_copy=X_I_indexes.copy()
      H=Gillespie_Direct_Method_Non_markovian_network(N=N,X_I_indexes=X_I_indexes_copy,method='aproximado',params=[],densities=[psi_infection,psi_recovery],survivals=[psi_survival_infection,psi_survival_recovery],G_network_initial=G_copy,T=14)
    if(k==0):
        plt.plot(H['time'],H['I(t)'], 'orange',linewidth=0.7,label='General aproximado (MG-T)')
    else:
        plt.plot(H['time'],H['I(t)'], 'orange',linewidth=0.7)
    i=np.argmax(H['I(t)'])
    picos['genApprox'].append([H['time'][i],H['I(t)'][i]])
    times_general_approx.append((time.time() - iter_time))
    matriz_trayectorias.append(H['time'])
    matriz_trayectorias.append(H['I(t)'])
    matriz_trayectorias.append(['MW'])
    #------------------------------------
    #------------------------------------ Weibull
    iter_time = time.time()
    H=[]
    while(len(H)<100):
      G_copy=G.copy()
      X_I_indexes_copy=X_I_indexes.copy()
      H=Gillespie_Direct_Method_Non_markovian_network(N=N,X_I_indexes=X_I_indexes_copy,method='Weibull',params=[a_I,a_R,lammbda_I,lammbda_R],G_network_initial=G_copy,T=14)
    if(k==0):
        plt.plot(H['time'],H['I(t)'], 'purple',linewidth=0.7,label='Weibull (MW)')
    else:
        plt.plot(H['time'],H['I(t)'], 'purple',linewidth=0.7)
    i=np.argmax(H['I(t)'])
    picos['weib'].append([H['time'][i],H['I(t)'][i]])
    times_weibull.append((time.time() - iter_time))
    matriz_trayectorias.append(H['time'])
    matriz_trayectorias.append(H['I(t)'])
    matriz_trayectorias.append(['MM'])
    #------------------------------------
    #------------------------------------ Markoviano
    iter_time = time.time()
    H=[]
    while(len(H)<100):
      G_copy=G.copy()
      X_I_indexes_copy=X_I_indexes.copy()
      H=Gillespie_Direct_Method_Network(beta,gamma,N,X_I_indexes=X_I_indexes_copy,G_network_initial=G_copy,T=14)
    if(k==0):
        plt.plot(H['time'],H['I(t)'], 'black',linewidth=0.7,label='Markoviano (MM) ')
    else:
        plt.plot(H['time'],H['I(t)'], 'black',linewidth=0.7)
    i=np.argmax(H['I(t)'])
    picos['mark'].append([H['time'][i],H['I(t)'][i]])
    times_markovian.append((time.time() - iter_time))
    matriz_trayectorias.append(H['time'])
    matriz_trayectorias.append(H['I(t)'])
    matriz_trayectorias.append(['MG'])
    #------------------------------------
plt.legend()
plt.savefig('comparacion.png',dpi=300)
plt.show()
print("--- %s Segundos en promedio para ejecutar 1 iteracion del general ---" %np.mean(times_general))
print("--- %s Segundos en promedio para ejecutar 1 iteracion del general usando aproximacion---" %np.mean(times_general_approx))
print("--- %s Segundos en promedio para ejecutar 1 iteracion del caso Weibull ---" % np.mean(times_weibull))
print("--- %s Segundos en promedio para ejecutar 1 iteracion del caso markoviano ---" % np.mean(times_markovian))
