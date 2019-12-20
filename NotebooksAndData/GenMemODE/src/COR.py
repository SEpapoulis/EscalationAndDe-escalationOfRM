from scipy.stats import uniform
import numpy as np
import pandas as pd
import multiprocessing
import os

from GenMemODE import General_interaction,Parallel_interaction,Memory_interaction,load_traits,Sr_sweep,Producer_Consumer_ODE,GenMem_ODE,find_endpoints


def pool_run(results,Sr_rng,model,initvalues,traits,t_final,steps,seed):
    dfs = Sr_sweep(Sr_rng,model,initvalues,traits,t_final,steps,seed,print_status=False)
    df=find_endpoints(dfs)
    df.to_csv(results,index=False)
    
def get_init(members):
    initvalues = [('R',1)]
    for member in members:
        initvalues.append(tuple([member,100]))
    return(np.array(initvalues))


def COR_Analysis(parameters,results_dir):
    args=[]
    orgRM=['','A']
    Sr_rng = np.logspace(4.5,8,100)

    t_final = 2000
    steps = 500

    for index,row in parameters.iterrows():
        i=str(index)

        #General
        f = os.path.join(results_dir,"Gen_rep_{}.csv".format(i))
        r,members,c=General_interaction(orgRM,
                                        resistance_function=lambda x:(row['r']**x),
                                        cost_function=lambda x:(row['c'])**x)
        traits = load_traits(row['alpha'],row['phi'],row['delta_p'],row['delta_c'],row['beta'],c,r)
        args.append([f,Sr_rng,Producer_Consumer_ODE,get_init(members),traits,t_final,steps,True])

        #Parallel
        f = os.path.join(results_dir,"Par_rep_{}.csv".format(i))
        r,members,c=Parallel_interaction(orgRM,
                                        resistance_function=lambda x:(row['r']**x),
                                        cost_function=lambda x:(row['c'])**x)
        traits = load_traits(row['alpha'],row['phi'],row['delta_p'],row['delta_c'],row['beta'],c,r)
        args.append([f,Sr_rng,Producer_Consumer_ODE,get_init(members),traits,t_final,steps,True])

        #Memory
        f = os.path.join(results_dir,"Mem_rep_{}.csv".format(i))
        r,members,c=Memory_interaction(orgRM,
                                        resistance_function=lambda x:(row['r']**x),
                                        cost_function=lambda x:(row['c'])**x,
                                        partial_resistance_function=lambda pr:(.9**pr))
        traits = load_traits(row['alpha'],row['phi'],row['delta_p'],row['delta_c'],row['beta'],c,r)
        args.append([f,Sr_rng,GenMem_ODE,get_init(members),traits,t_final,steps,True])
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()-1) as pool:
        pool.starmap(pool_run,args)

m=[]
alpha,phi,delta_p,delta_c,beta=1.0,1e-8,.2,.2,25.0
for r in [.001,.01,.5]:
    for c in np.linspace(0.01,.99,99):
        m.append([alpha,phi,delta_p,delta_c,beta,c,r])
parameters = pd.DataFrame(m)
parameters.to_csv('data/COR_parameters.csv')
parameters.columns = ['alpha','phi','delta_p','delta_c','beta','c','r']
COR_Analysis(parameters,'data/COR/')