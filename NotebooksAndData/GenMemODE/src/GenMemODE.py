import numpy as np
import pandas as pd
from scipy import integrate

def GenMem_ODE(state,t,Sr,traits):
    r'''
    General to Memory interaction Model

    For this model to work, we must ALWAYS assume the r matrix is ixj, where i == j
    in this model, beta CANNOT vary betwen viruses

    .. math::
        \frac{dR}{dt} = Sr - \sum_{i=1}^{n} \alpha_{i} (1-c_{i}) P_{i} R
    .. math::
        \frac{dP_{i}}{dt} = \alpha_{i} (1-c_{i}) P_{i} R -  \sum_{j=1}^{n} \phi_{i} r_{ji} P_{i} C_{j} - \delta_{p} P_{i}
    .. math::
        \frac{dC_{j}}{dt} = \beta_{j} \sum_{i=1}^{n} \phi_{i} r_{ji} P_{i} C_{j} - \delta_{c} C_{j}

    '''

    R,orgs = state[0],state[1:]

    p_index = np.where(traits['is_p']==1)[0]

    c_index = np.where(traits['is_p']==0)[0]
    
    pprod = calc_pprod(traits,R,p_index,orgs)
    ppred = calc_ppred(traits,p_index,c_index,orgs)
    ploss = calc_ploss(traits,p_index,orgs)
    cprod = calc_cprod_mem(traits,p_index,c_index,orgs) #mass balance with ppred
    closs = calc_closs(traits,c_index,orgs)
    
    #ODE
    #TODO: implement memory at this stage, where memory can shift dCdt[1:] to dCdt[0]
    dRdt = Sr - np.sum(pprod)
    dPdt = pprod - ppred - ploss
    dCdt = cprod - closs
    return np.concatenate([np.r_[[dRdt]],dPdt,dCdt])

#TODO
def calc_cprod_mem(traits,p_index,c_index,orgs):
    r'''
    Consumer Production, with memory

    .. math::
        \beta_{j} \sum_{i=1}^{n} \phi_{i} r_{ji} P_{i} C_{j} 
    '''
    producers = orgs[p_index]
    consumers = orgs[c_index]
    phis=traits['phi']
    beta=traits['beta']
    #Because production of each virus can only come from a single producer, we must do our matrix
    #with producers as columns
    return(phis*producers*np.dot(consumers*beta,traits['r'].T))


def Producer_Consumer_ODE(state,t,Sr,traits):
    r'''
    Producer Consumer Interaction Model

    .. math::
        \frac{dR}{dt} = Sr - \sum_{i=1}^{n} \alpha_{i} (1-c_{i}) P_{i} R
    .. math::
        \frac{dP_{i}}{dt} = \alpha_{i} (1-c_{i}) P_{i} R -  \sum_{j=1}^{n} \phi_{i} r_{ji} P_{i} C_{j} - \delta_{p} P_{i}
    .. math::
        \frac{dC_{j}}{dt} = \beta_{j} \sum_{i=1}^{n} \phi_{i} r_{ji} P_{i} C_{j} - \delta_{c} C_{j}

    '''

    R,orgs = state[0],state[1:]

    p_index = np.where(traits['is_p']==1)[0]

    c_index = np.where(traits['is_p']==0)[0]
    
    pprod = calc_pprod(traits,R,p_index,orgs)
    ppred = calc_ppred(traits,p_index,c_index,orgs)
    ploss = calc_ploss(traits,p_index,orgs)
    cprod = calc_cprod(traits,p_index,c_index,orgs) #mass balance with ppred
    closs = calc_closs(traits,c_index,orgs)
    
    #ODE
    dRdt = Sr - np.sum(pprod)
    dPdt = pprod - ppred - ploss
    dCdt = cprod - closs
    return np.concatenate([np.r_[[dRdt]],dPdt,dCdt])

def calc_pprod(traits,R,p_index,orgs):
    r'''
    Producer Production

    .. math::
       \alpha_{i} (1-c_{i}) P_{i} R
    '''
    
    alphas = traits['alpha'][p_index]
    costs = traits['c'][p_index]
    producers = orgs[p_index]#builds a new array to be used specifically for the following computation	
    prod = alphas*(1-costs)*producers*R
    return(prod)

def calc_ppred(traits,p_index,c_index,orgs):
    r'''
    Producer Predation

    .. math::
        \sum_{j=1}^{n} \phi_{i} r_{ji} P_{i} C_{j}
    '''
    producers = orgs[p_index]
    consumers = orgs[c_index]
    phis=traits['phi']
    #we need to transpose the r matrix from rows of P and columns of C to rows of C and columns of P
    #to get the correct dot product
    return(producers*phis*np.dot(consumers,traits['r'].T))
    
def calc_ploss(traits,p_index,orgs):
    r'''
    Producer Loss

    .. math::
        \delta_{p} P_{i}
    '''
    delta_p = traits['delta_p'][p_index]
    producers = orgs[p_index]
    return(producers*delta_p)

def calc_cprod(traits,p_index,c_index,orgs):
    r'''
    Consumer Production

    .. math::
        \beta_{j} \sum_{i=1}^{n} \phi_{i} r_{ji} P_{i} C_{j} 
    '''
    producers = orgs[p_index]
    consumers = orgs[c_index]
    phis=traits['phi']
    beta=traits['beta']
    #we need to transpose the r matrix from rows of P and columns of C to rows of C and columns of P
    #to get the correct dot product
    return(beta*consumers*np.dot(producers*phis,traits['r']))

def calc_closs(traits,c_index,orgs):
    r'''
    Consumer Loss

    .. math::
        \delta_{c} C_{j} 
    '''
    delta_c = traits['delta_c']#[c_index]
    consumers = orgs[c_index]
    return(consumers*delta_c)

def General_interaction(host_defense,resistance_function=lambda r:(.01**r),cost_function=lambda c:(.1*c)):
    r'''
    General consumer producer interaction

    This consumer/producer interaction assumes one consumer can consume all producers. 
    If using hosts and viruses, producers = hosts while consumers = viruses

    Parameters
    ----------
    host_defense : list of str
        a list of strings, where each character in the string represents a defense system.
        This model assumes that the resistance confered by each defense system is identical.
    resistance_function : function
        resistance function used to calculate total resistance of each host. Default: r = .001^(num defense systems)

    Returns
    -------
    numpy array
        rows are the producers while columns are the consumers
    members
        a list of the members that will be in the simulation (including the general consumer at index -1)
    numpy array
        a numpy array of costs
    '''
    interaction=[]
    members=[]
    costs=[]
    for i,RM in enumerate(host_defense):
        interaction.append([resistance_function(len(set(RM)))])
        members.append('P_{}'.format(str(i)))
        costs.append(cost_function(len(set(RM))))
    members.append('C')
    return(np.array(interaction),members,np.array(costs))

def Parallel_interaction(host_defense,resistance_function=lambda r:(.01**r),cost_function=lambda c:(.1*c)):
    r'''
    Parallel consumer producer interaction

    This consumer/producer interaction assumes one consumer for each producer. 
    If using hosts and viruses, producers = hosts while consumers = viruses

    Parameters
    ----------
    host_defense : list of str
        a list of strings, where each character in the string represents a defense system.
        This model assumes that the resistance confered by each defense system is identical.
    resistance_function : function
        resistance function used to calculate total resistance of each host. Default: r = .001^(num defense systems)
    cost_function : function
        cost function used to calculate the cost of resistance. Default: c = 0.1*(num defense sysetms)
    Returns
    -------
    numpy array
        rows are the producers while columns are the consumers
    members
        a list of the members that will be in the simulation (including the general consumer at index -1)
    numpy array
        a numpy array of costs
    '''
    interaction=[]
    members_p=[]
    members_c=[]
    costs=[]
    for i,RM in enumerate(host_defense):
        row = [0] *len(host_defense)
        costs.append(cost_function(len(set(RM))))
        row[i] = resistance_function(len(set(RM)))
        interaction.append(row)
        members_p.append('P_{}'.format(str(i)))
        members_c.append('C_{}'.format(str(i)))
    return(np.array(interaction),members_p+members_c,np.array(costs))

def Memory_interaction(host_defense,resistance_function=lambda r:(.01**r),cost_function=lambda c:(.1*c),partial_resistance_function=lambda pr:(1**pr)):
    r'''
    Memory consumer producer interaction

    This consumer/producer interaction assumes that consumers can "remember" the previous
    host that was infected, thus, the identity of the host defense system is critical in 
    determining reistance. The total number of "effective defense sysetms" will be used
    to calculate resistance from the resistance_function. "effecitve defense sysetms" are
    defined as difference in the host defense sysetms and the defense systems of the previous
    host infected. For example, host 'ABC' and phage 'ABC' results in zero conferred resistance 
    while host 'ABD' would have the effective defense of 1 system

    Parameters
    ----------
    host_defense : list of str
        a list of strings, where each character in the string represents a defense system.
        This model assumes that the resistance confered by each defense system is identical.
    resistance_function : function
        resistance function used to calculate total resistance of each host. Default: r = .001^(num defense systems)
    cost_function : function
        cost function used to calculate the cost of resistance. Default: c = 0.1*(num defense sysetms)

    Returns
    -------
    numpy array
        rows are the producers while columns are the consumers
    members
        a list of the members that will be in the simulation (including the general consumer at index -1)
    numpy array
        a numpy array of costs
    '''
    interaction=[]
    members_p=[]
    members_c=[]
    costs = []
    for i,RM_producer in enumerate(host_defense):
        row=[]
        costs.append(cost_function(len(RM_producer)))
        for RM_consumer in host_defense:
            effective_resistance = len(set(RM_producer) - set(RM_consumer))
            partial_resistance = len(set(RM_producer).intersection(RM_consumer))
            resistance = resistance_function(effective_resistance) * partial_resistance_function(partial_resistance)
            row.append(resistance)
        interaction.append(row)
        members_p.append('P_{}'.format(str(i)))
        members_c.append('C_{}'.format(str(i)))
    return(np.array(interaction),members_p+members_c,np.array(costs))

def load_traits(alpha,phi,delta_p,delta_c,beta,c,r_interaction,additional = []):
    r'''
    Loads the traits needed to run simulations

    if floats are passed through alpha, phi, delta_p,delta_c, beta, and c, all producers will
    be assigned that parameter value. In order to specify a paramemter value, list/arrays of 
    proper length must be passed. "i" is the number of producers while "j" is the number of consumers

    Parameters
    ----------
    alpha : float or array of length i
        Resource utilization
    phi : float or array of length i
        Baseline resistance of producers to consumers
    delta_p : float or array of length i
        Loss of producers
    delta_c : float or array of length j
        Loss of consumers
    beta : float or array of length j
        Consumer produced per producer preyed upon. Commonly known as "Burst size" when consumers are viruses
    c : float or array of length i
        Cost of a producers resistance
    r_interaction : 2Darray
        rows are producers, columns are consumers. element [i,j] would be the resistance of producer i to consumer j
    additional : list of tuples, optional
        list of tuples, where tuples consist of the trait name and numpy array of values

    Return
    ------
    Dictionary of Traits

    '''

    num_producers = len(r_interaction[:,0])
    num_consumers = len(r_interaction[0])
    if type(alpha) == float or type(alpha) == np.float64:
        alpha = np.full(num_producers,alpha)
    if type(phi) == float or type(phi) == np.float64:
        phi = np.full(num_producers,phi)
    if type(delta_p) == float or type(delta_p) == np.float64:
        delta_p = np.full(num_producers,delta_p)
    if type(delta_c) == float or type(delta_c) == np.float64:
        delta_c = np.full(num_consumers,delta_c)
    if type(beta) == float or type(beta) == np.float64:
        beta = np.full(num_consumers,beta)
    if type(c) == float or type(c) == np.float64:
        c = np.full(num_producers,c)

    traits = {'alpha' : alpha,
                'phi' : phi,
                'delta_p' : delta_p,
                'delta_c' : delta_c,
                'beta' : beta,
                'c' : c,
                'r' : r_interaction,
                'is_p' : ( np.concatenate([np.ones(num_producers,dtype=int),np.zeros(num_consumers,dtype=int)]) )}
    if additional:
        for name,array in additional:
            traits[name] = array
    
    return(traits)


def Launch_Numerical(model,Sr,traits,init,t_end,steps_per_t=1):
    t_series=np.linspace(0,t_end,t_end*steps_per_t)
    init_values=[]
    result_columns=['t']
    for member,val in init:
        init_values.append(val)
        result_columns.append(member)
    #Integration
    raw = integrate.odeint(model, y0=init_values, t=t_series,args=(Sr,traits))
    results = pd.DataFrame(raw).reset_index()
    results.columns = result_columns
    #Adding timesteps per t
    results['t'] = results['t'] / steps_per_t
    #removing timesteps that are not whole numbers
    results = results[results['t']%1 == 0]
    return(results)

def Sr_sweep(Srs,model,initvalues,traits,t_final,steps,seed=True,print_status=True):
    '''
    Iterativly launches numerical simulations with different Srs
    
    Parameters
    ----------
    Srs : array
        An array indicating the Sr of each simulation
    model : function
        Model used in numerical integration
    initvalues : list of tuples
        a list of tuples, where tuple[0]= member name and tuple[1]= inital value in simulation
    traits: dict of arrays
        A dictionary mapping trait names to arrays. Note that this must be compatible with the respective model
    t_final : int
        How long the simulation should run
    steps : int
        How many steps should be taken per t
    seed : bool, optional
        Seed the initial values of the next simulation with the results of the previous
    
    '''
    
    results = []
    init_col = [el[0] for el in initvalues]
    num_sim = len(Srs)
    if print_status:
        print("Preparing to run {} simulations between {} and {}".format(str(num_sim),min(Srs),max(Srs)))
    for i,Sr in enumerate(Srs):
        if print_status:
            print("{:.2f}% Complete".format(i/num_sim*100),end='\r')
        temp = Launch_Numerical(model=model,Sr=Sr,traits=traits,init=initvalues,t_end=t_final,steps_per_t=steps)
        temp['Sr']=Sr
        if seed:
            last=np.array(temp.iloc[-1][init_col])
            #rebuild as list of tuples
            initvalues=[]
            for i,col in enumerate(init_col):
                if last[i] < 1:
                    last[i]=1
                initvalues.append((col,last[i]))
        results.append(temp)
    if print_status:
        print("100.00% Complete")
    return(results)

def find_endpoints(results):
    '''
    Build a dataframe composed of the last index from a list of dataframes
    '''
    m=[]
    for el in results:
        m.append(el.iloc[-1])
    return(pd.DataFrame(m))