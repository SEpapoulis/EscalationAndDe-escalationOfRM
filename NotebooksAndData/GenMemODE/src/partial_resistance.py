from scipy.stats import uniform
import numpy as np
import pandas as pd
import multiprocessing
import os

from GenMemODE import General_interaction,Parallel_interaction,Memory_interaction,load_traits,Sr_sweep,Producer_Consumer_ODE,GenMem_ODE,find_endpoints




    
def get_init(members):
    initvalues = [('R',1)]
    for member in members:
        initvalues.append(tuple([member,100]))
    return(np.array(initvalues))

def pool_run(results,Sr_rng,model,initvalues,traits,t_final,steps,seed):
    print('starting sim')
    dfs = Sr_sweep(Sr_rng,model,initvalues,traits,t_final,steps,seed,print_status=False)
    df=find_endpoints(dfs)
    df.to_csv(results,index=False)


def partial_resistance(parameters,orgRM,results_dir):
    args=[]
    Sr_rng = np.logspace(4.5,8,100)

    t_final = 2000
    steps = 500
    for index,row in parameters.iterrows():
        i=str(index)
        #Memory
        f = os.path.join(results_dir,"Mem_rep_{}.csv".format(i))
        r,members,c=Memory_interaction(orgRM,
                                        partial_resistance_function=lambda pr:(row['pr']**pr))
        traits = load_traits(row['alpha'],row['phi'],row['delta_p'],row['delta_c'],row['beta'],c,r)
        args.append([f,Sr_rng,GenMem_ODE,get_init(members),traits,t_final,steps,True])
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()-1) as pool:
        pool.starmap(pool_run,args)

alpha,phi,delta_p,delta_c,beta=1.0,1e-8,.2,.2,25.0
prs = np.around(np.arange(start=0.01,stop=.51,step=.01),4)

parameters = {'pr':prs,
            'alpha':np.array([alpha]*len(prs)),
            'phi':np.array([phi]*len(prs)),
            'delta_p':np.array([delta_p]*len(prs)),
            'delta_c':np.array([delta_c]*len(prs)),
            'beta':np.array([beta]*len(prs)),
            'c':np.array([.1]*len(prs))

}
parameters = pd.DataFrame(parameters)

orgRM=[]
rmcounter=0
for i in range(0,10):
    orgRM.append([])
    while i > len(orgRM[-1]):
        orgRM[-1].append(rmcounter)
        rmcounter+=1
print('Starting unique')
partial_resistance(parameters,orgRM,'data/partial_r/orgRM_unique')
print('Unique Complete')
orgRM=[]
rmcounter=0
for i in range(0,10):
    orgRM.append([])
    rmcounter=0
    while i > len(orgRM[-1]):
        orgRM[-1].append(rmcounter)
        rmcounter+=1
print('Starting stacking')
partial_resistance(parameters,orgRM,'data/partial_r/orgRM_stacking')
print('Stacking Complete')