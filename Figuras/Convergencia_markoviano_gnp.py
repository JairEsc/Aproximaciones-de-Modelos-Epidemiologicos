#Ejemplo
import time
colores=['green','orange','purple']
np.random.seed(0)
for pot in [2,3,4]:
    iter_time = time.time()
    beta=75
    gamma=1
    N=10**pot
    ksum=1
    G=nx.fast_gnp_random_graph(n=N, p=0.1, seed=None, directed=False)
    X_I_indexes=np.random.choice(N,size=int(N*0.01),replace=False).tolist()
    for k in range(N):
        if k in X_I_indexes:
            G.nodes[k]['Estado:']='Infeccioso'
        else:
            G.nodes[k]['Estado:']='Susceptible'
    print("Empieza")
    #------------------------------------ Markoviano
    for iter in range(5):
        H=[]
        while(len(H)<10):
          G_copy=G.copy()
          X_I_indexes_copy=X_I_indexes.copy()
          H=Gillespie_Direct_Method_Network(beta,gamma,N,X_I_indexes=X_I_indexes_copy,G_network_initial=G_copy,T=6)
        if(iter==0):
            plt.plot(H['time'],H['I(t)']/N, colores[pot-2],linewidth=0.4,label='N=10e'+str(pot))
        else:
            plt.plot(H['time'],H['I(t)']/N, colores[pot-2],linewidth=0.4)
        #------------------------------------
    print("Tomó "+str(time.time()-iter_time)+'s. para 5 iteraciones con N=10e'+str(pot))
INPUT = (0.99*N, 0.01*N, 0.0)
t_start = 0.0; t_end = 6; t_inc = .01
t_range = np.arange(t_start, t_end+t_inc, t_inc)
plt.plot(t_range,g(t_range,INPUT,[beta*0.1,gamma,N])/N, 'k',label='Campo Medio')
plt.legend()
plt.savefig('markoviano_convergencia_gnp.png',dpi=300)
plt.show()
