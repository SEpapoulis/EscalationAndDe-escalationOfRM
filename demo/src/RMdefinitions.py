MTaseDomains=set(['DNA_methylase','Eco57I','EcoRI_methylase','HsdM_N','MethyltransfD12','MT-A70',
'N6_Mtase','N6_N4_Mtase','TypeIII_RM_meth','Dam'])

#part of Type IV restriction enzymes. May indicate only digest modified DNA
MrrDomains=set(['Mrr_cat_2','Mrr_cat','Mrr_N'])
# Type I specificity subunit
Type_IS=set(['Methylase_S'])

#Domains found in Type I systems
Type_I=['HsdM_N','HSDR_N','HSDR_N_2','EcoR124_C','EcoEI_R_C']

#Domains found in Type III systmes
Type_III=['TypeIII_RM_meth','ResIII']

#Domains responsible for defining target recognition sequence of MTase/REase
recognition_domains = set(['Methylase_S','TaqI_C'])

#if a methyltransferase is this big, there is a chance there is a undetected endonuclease domain 
min_pTIIG_size=750 

#hypothesized to be endonuclease?
Hypothesized_Endonucleases=set(['Uma2'])

def process_exception(hit):
    if '.' in hit:
        if hit[0] == 'M':
            return('M')
        elif hit[0] == 'S':
            return('S')
        else:
            return('R')
    else:
        return('R')


###
#ALL DOMAINS NEED TO BE IN SETS
###
def is_hypendo(domains):
    if domains == Hypothesized_Endonucleases:
        return(True)
    else:
        return(False)

def has_TRD(domains):
    common = domains.intersection(recognition_domains)
    if common:
        return(True)
    return(False)

def has_MTase(domains):
    common = domains.intersection(MTaseDomains)
    if common:
        return(True)
    return(False)

def has_mT3(domains):
    common = domains.intersection(set(['TypeIII_RM_meth']))
    if common:
        return(True)
    return(False)

def has_REase(domains):#excludes Mrr
    noMT = domains.difference(MTaseDomains)
    noMT_noTRD = noMT.difference(recognition_domains)
    noMT_noTRD_noTIV = noMT_noTRD.difference(MrrDomains)
    if noMT_noTRD_noTIV:
        return(True)
    return(False)

def has_T1(domains):
    common = domains.intersection(Type_I)
    if common:
        return(True)
    return(False)

def has_TIV(domains):
    common = domains.intersection(MrrDomains)
    if common:
        return(True)
    return(False)

def only_TIS(domains):
    if domains == Type_IS:
        return(True)
    return(False)

def is_rT3(domains):
    if domains == set(['ResIII']):
        return(True)
    return(False)

def is_putativeTIIG(length):
    if length >= min_pTIIG_size:
        return(True)
    return(False)

def is_TIIG(domains,length):
    if is_putativeTIIG(length):
        if has_MTase(domains) and has_REase(domains): #has MTase and endo
            return(True)
        elif has_TIV(domains) and has_REase(domains): #has mrr with others
            return(True)
        else:
            return(False)
    else:
        return(False)
