from RMdefinitions import *
from datastructures import Fasta
import pandas as pd
import numpy as np

def load_simpletxt(file):
    dat=[]
    with open(file,'r') as f:
        for line in f:
            dat.append(line.strip())
    return(dat)

#Change if directory names are altered!###########################################
resultsdir = 'data/RMsearch/results/'
RMprots=Fasta('data/RMsearch/results/ElementsFound.fasta')
noRMorgs = set(load_simpletxt("data/RMsearch/results/No_element_orgs.txt"))
failed = set(load_simpletxt("data/RMsearch/results/failed.txt"))
genome_stats = pd.read_csv("data/RMsearch/results/orgid_dat.csv",header=None)
genome_stats.columns = ['assembly', 'GenomeSize','NumContigs']
assemblydat = pd.read_csv("data/RMsearch/demo_assemblies.csv")
##########################################################


def classify(dat):
    domains,length,product,blastsummary=dat
    domains = set(domains)
    length=int(length)
    typeclass=str()
    if (not domains) and not pd.isnull(blastsummary):
        code=process_exception(blastsummary[0])
        if code == 'R':
            typeclass='rT2'
        elif code == 'S':
            typeclass='sT1'
        else:
            typeclass='mT2'

    elif has_MTase(domains):
        if is_TIIG(domains,length):
            typeclass='T2G'
        elif has_T1(domains):
            typeclass='mT1'
        elif has_mT3(domains):
            typeclass='mT3'
        elif is_putativeTIIG(length):
            typeclass='pT2G'
        else:
            typeclass='mT2'
    
    elif only_TIS(domains):
        typeclass='sT1'

    #Checking for case where there is a detectable Mrr_cat domain and
    #other Rease domains
    elif is_TIIG(domains,length):
        typeclass='T2G'
        
    elif has_TIV(domains):
        typeclass='T4'

    elif has_REase(domains):
        if has_T1(domains):
            typeclass='rT1'
        elif is_rT3(domains):
            typeclass='rT3'
        elif is_hypendo(domains):
            typeclass='pr'
        else:
            typeclass='rT2'
    
    #if it recognizes DNA and does not have MTase domain, must be endonuclease?
    elif has_TRD(domains):
        typeclass = 'rT2'

    else:
        print('WARNING, NEVER CLASSIFIED')
        print('\t'.join([str(domains),length,product,str(blastsummary)]))
        typeclass=''
    
    return(typeclass)

def get_RMcount(RMtypedict,RMtype=str()):
    m = RMtypedict['mT'+RMtype]
    r = RMtypedict['rT'+RMtype]
    if m > r:
        return(r)
    else:
        return(m)
    

#checks if accession is not a pseudo (< or > means truncation, see refseq documentation)
def not_pseudo(accessions):
    out=[]
    for accession in accessions:
        if '>' in accession or '<' in accession:
            out.append(False)
        else:
            out.append(True)
    return(out)

#removes domains from dataframe if all domains are to be ignored
def filter_domain(df,ignore=set()):
    out=[]
    for row in df.itertuples(index=False):
        domain = row[1]
        blastexcept = row[-1]
        #if there is not a domain, entry is the result of blasthit
        if pd.isnull(domain):
            out.append(True)
        #There is a domain
        else:
            for el in domain.split(';'):
                keep=False
                #make sure that all domains can be ignored and there is not a blastexcept
                #if there is a blastexcept, this statement is always true
                if el not in ignore or not pd.isnull(blastexcept):
                    keep=True
            out.append(keep)
    return(out)

#helper function for accession_vicinity
def in_vicinity(se1,se2,distance=4000):
    #is gene se1 within 4000 bp or se2
    s1,e1 = se1
    s2,e2 = se2

    if abs(s1 - e2 ) <= distance or abs(s2 - e1) <=distance:
        return(True)
    else:
        return(False)

#helper function for accession_vicinity
def _get_location(locstr):
    #returns a list of tupples. if there is no join, list length is 1
    loc_list=[]
    outstr=str()
    for char in locstr:
        if char.isdigit() or char == '.' or char ==',':
            outstr=outstr+char
    str_list = outstr.split(',')
    for el in str_list:
        _ = el.split('..')
        if len(_) == 2:
            loc_list.append( (int(_[0]) , int(_[1]) ) )
        else: #for when there is just one base pair that needs to be joined
            loc_list.append( (int(_[0]) , int(_[0]) ) )
    return(loc_list)

def accession_vicinity(accession1, accession2, dist = 4000):
    locus1,location1=accession1.split('|')
    locus2,location2=accession2.split('|')
    loc1=_get_location(location1)
    loc2=_get_location(location2)

    if locus1 != locus2: #not on the same DNA mol, not in vicinity
        return(False)

    if len(loc1) == 1 and len(loc2) == 1: #genes are not joined together
        if in_vicinity(loc1[0],loc2[0],dist):
            return(True)

    elif len(loc1) > 1 and len(loc2) == 1: #loc1 has join
        for el in loc1:
            if in_vicinity(el,loc2[0],dist):
                return(True)

    elif len(loc1) == 1 and len(loc2) > 1: #loc2 has join
        for el in loc2:
            if in_vicinity(el,loc1[0],dist):
                return(True)

    elif len(loc1) > 1 and len(loc2) > 1: #loc1 and loc2 have join, compare all
        for el1 in loc1:
            for el2 in loc2:
                if in_vicinity(el1,el2,dist):
                    return(True)
    return(False)

def localize_clusters(ordered_accesionlist,distance=4000,anchor_accessions=set()):
    #find genes that are localized 'distance' apart and splits them into sublists
    #where each list is a localized set of detected genes
    listlen=len(ordered_accesionlist)#length of RM results file
    if listlen <= 1:
        return([ordered_accesionlist])
    i=1
    out=[]
    temp=[ordered_accesionlist[0%listlen]]
    while i < listlen:
        accession=ordered_accesionlist[i%listlen]
        i_temp=i-1
        last_accession=ordered_accesionlist[i_temp%listlen]
        if accession_vicinity(accession,last_accession,distance):
            temp.append(accession)
        else:
            out.append(temp)
            temp=[accession]
        i+=1
    else:
        out.append(temp)
        #check to see if we need to merge index 0 and -1 (incase of circular Chromosome/plasmids)
        if len(out) >1 and accession_vicinity(accession,out[0][0],distance):
            temp=out.pop(0)
            out[-1] = out[-1] + temp


    if anchor_accessions: #looking for clusters with particular accessions
        newout=[]
        for sublist in out:
            keep = False
            for accession in sublist:
                if accession in anchor_accessions:
                    keep=True
            if keep:
                newout.append(sublist)
        out=newout
    return(out)


#building helper dictionaries
assembly_genomesize =  pd.Series(genome_stats.GenomeSize.values,index=genome_stats.assembly).to_dict()
assembly_NumContigs =  pd.Series(genome_stats.NumContigs.values,index=genome_stats.assembly).to_dict()
assembly_genus =  dict()
assembly_phylum =  dict()
for row in assemblydat.itertuples(index=False):
    assembly = row[0].split('.')[0]
    assembly_genus[assembly]=row[3]
    assembly_phylum[assembly]=row[2]

#recording RMstats for each org
org_RMdat=[]
#recording RM classificaiton for each protein found
RMprotclass=[]
#Fastsa file of all pTIIG RM systems
#Note pTIIG == protein with alignment to methyltransferase domain that is >750AA
pTIIG = Fasta()


org_RMdat={'orgid':[],'mT1':[],'rT1':[],'sT1':[],'rmT1':[],
            'mT2':[],'rT2':[],'rmT2':[],
            'mT3':[],'rT3':[],'rmT3':[],
            'T4':[],'T2G':[],'pT2G':[],
            'pr':[],'prm':[],'rmT?':[]}

#settings#####################
distance=4000
ignore_pseudo=True
ignore_domains = set(['Uma2'])
##############################
itot = len(assemblydat['assembly'])

#adding data for organisms without RM proteins
for assembly in noRMorgs:
    for el in org_RMdat:
        if el == 'orgid':
            org_RMdat[el].append(assembly)
        else:
            org_RMdat[el].append(0)

for i,assembly in enumerate(assemblydat['assembly']):
    pcnt = str(round((i+1)/itot*100,3))
    gcf = assembly.split('.')[0]
    print("Progress: {}%    Processing assembly: {}".format(pcnt,gcf),end='\r')
    #skip if there is noRMorgs, NOTE: GOES TO TOP OF FOR LOOP FOR NEXT ASSEMBLY
    if gcf in noRMorgs:
        for el in org_RMdat:
            if el == 'orgid':
                org_RMdat[el].append(gcf)
            else:
                org_RMdat[el].append(0)
        continue

    #skiping files that failed, NOTE: GOES TO TOP OF FOR LOOP FOR NEXT ASSEMBLY
    elif gcf in failed:
        continue

    org_results = pd.read_csv(resultsdir+gcf+'.bioscan.csv')

    if ignore_pseudo:
        org_results = org_results[not_pseudo(org_results['accession'])]

    if ignore_domains:
        org_results=org_results[filter_domain(org_results,ignore_domains)]

    #Error check to make sure org_results is not empty, NOTE: GOES TO TOP OF FOR LOOP FOR NEXT ASSEMBLY
    if org_results.empty:
        for el in org_RMdat:
            if el == 'orgid':
                org_RMdat[el].append(gcf)
            else:
                org_RMdat[el].append(0)
        continue
    #getting assessions that are MTase positive
    MtaseAccessions=list()
    for row in org_results.itertuples(index=False):
        #checking domain
        if row[1] in MTaseDomains:
            MtaseAccessions.append(row[0])
        #checking BlastException
        elif not pd.isnull(row[-1]):
            if row[-1][0:2] == 'M.':
                MtaseAccessions.append(row[0])
                    
    #Generating clusters of hits next to each other in the chromosome
    clusters = localize_clusters(list(org_results['accession']),distance,anchor_accessions=MtaseAccessions)
    cluster_accessions={}
    clustered=list()
    for i,cluster in enumerate(clusters):
        temp=[]
        if len(cluster) >1:
            cluster_accessions[i] = cluster
            clustered.extend(cluster)
    clustered=set(clustered)
    not_clustered = set(org_results['accession'])-clustered
    org_results=org_results.set_index('accession')

    #counter for current organism
    RMtype={'mT1':0,'rT1':0,'sT1':0,'rmT1':0,
                'mT2':0,'rT2':0,'rmT2':0,
                'mT3':0,'rT3':0,'rmT3':0,
                'T4':0,'T2G':0,'pT2G':0,
                'pr':0,'prm':0,'rmT?':0}

    #evaluating proteins that were not in clusters

    for accession in not_clustered:
        dat=list(org_results.loc[accession])
        length=dat[1]
        if pd.isnull(dat[0]):
            dat[0]=''
        dat[0] = dat[0].split(';')

        prot_type = classify(dat)
        RMtype[prot_type]=RMtype[prot_type]+1
        ## adding TIIG to fasta for further annotation
        if prot_type == 'pT2G':
            pTIIG['>'+accession+'\t'+gcf]=RMprots[accession]
        #Adding to master list of identified proteins#########
        RMprotclass.append([gcf,accession,length,prot_type,np.nan,0,0,0,0])
            
    #We now need to evalute clusters of MTases and other proteins to see if
    #is an acceptable RM pair
    for cluster in cluster_accessions:
        temp_RMprotclass=[]
        cluster_RMtype={'mT1':0,'rT1':0,'sT1':0,'mT2':0,'rT2':0,'mT3':0,'rT3':0,
            'T4':0,'T2G':0,'pT2G':0,'pr':0}
        for accession in cluster_accessions[cluster]:
            dat=list(org_results.loc[accession])
            length=dat[1]
            if pd.isnull(dat[0]):
                dat[0]=''
            dat[0] = dat[0].split(';')
            
            prot_type = classify(dat)
            cluster_RMtype[prot_type]=cluster_RMtype[prot_type]+1
            ## adding TIIG to fasta for further annotation
            if prot_type == 'pT2G':
                pTIIG['>'+accession+'\t'+gcf]=RMprots[accession]

            #Adding to temporary master list of identified proteins#########
            temp_RMprotclass.append([gcf,accession,length,prot_type,cluster])
            
        #We will now alter II/III classifications based on context
        #If there is a type III endonuclease next to a methyltransferase,
        #we will alter the methyltransferase and call it a TypeIII methyltransferase
        #since not all recorded TypeIII RM systems are TypeIII_RM_meth domains
        while (cluster_RMtype['rT3'] > cluster_RMtype['mT3']) and cluster_RMtype['mT2']:
                cluster_RMtype['mT2'] = cluster_RMtype['mT2'] -1
                cluster_RMtype['mT3'] = cluster_RMtype['mT3'] +1
        #find the number of complete RM systems we can catagorize
        rmT1= get_RMcount(cluster_RMtype,'1')
        rmT2= get_RMcount(cluster_RMtype,'2')
        rmT3= get_RMcount(cluster_RMtype,'3')


        #find the number of orphans
        orphanR = cluster_RMtype['rT1'] + cluster_RMtype['rT2'] + cluster_RMtype['rT3'] - rmT3 - rmT2 - rmT1
        orphanM = cluster_RMtype['mT1'] + cluster_RMtype['mT2'] + cluster_RMtype['mT3'] - rmT3 - rmT2 - rmT1
        
        #check if orphans may be prm
        #note, prm would be rm systems classified with putative endonuclease
        #in the RMdefinitions file
        if orphanM >= cluster_RMtype['pr']:
            prm = cluster_RMtype['pr']
            orphanM = orphanM - prm
        else:
            prm = orphanM
            orphanM = orphanM - prm
            
        #if there are still r or m left, lets add them to an unclassified type
        if orphanM >=orphanR:
            rmUNK = orphanR
        else:
            rmUNK = orphanM

        
            
        for line in temp_RMprotclass:
            RMprotclass.append(line+[rmT1,rmT2,rmT3,prm])
            
        
        
        #adding cluster data to master dictionaries
        RMtype['rmT1']=RMtype['rmT1']+rmT1
        RMtype['rmT2']=RMtype['rmT2']+rmT2
        RMtype['rmT3']=RMtype['rmT3']+rmT3
        RMtype['prm']=RMtype['prm']+prm
        RMtype['rmT?']=RMtype['rmT?']+rmUNK
        for el in cluster_RMtype:
            RMtype[el]=RMtype[el]+cluster_RMtype[el]
            
    org_RMdat['orgid'].append(gcf)
    for key in RMtype:
        org_RMdat[key].append(RMtype[key])

phyla=[]
genus=[]
numcontigs=[]
genomesize=[]
for assembly in org_RMdat['orgid']:
    assembly=assembly.split('.')[0]
    try:
        genomesize.append(assembly_genomesize[assembly])
    except KeyError:
        genomesize.append('')
    try:
        numcontigs.append(assembly_NumContigs[assembly])
    except KeyError:
        numcontigs.append('')
    genus.append(assembly_genus[assembly])
    phyla.append(assembly_phylum[assembly])
#adding data to a master dictionary
org_RMdat['phyla'] = phyla
org_RMdat['genus'] = genus
org_RMdat['GenomeSize(bp)'] = genomesize
org_RMdat['NumContigs'] = numcontigs    
df = pd.DataFrame(org_RMdat)
df=df[['orgid','genus','phyla','GenomeSize(bp)','NumContigs','mT1','rT1','sT1','rmT1',
            'mT2','rT2','rmT2',
            'mT3','rT3','rmT3',
            'T4','T2G',
            'pr','prm','rmT?','pT2G']]    
#Saving Data
df.to_csv('data/RMsearch/RMdat.csv',index=False)
pTIIG.write_file('data/RMsearch/pTIIG/pTIIG.fasta')
RMprotclass=pd.DataFrame(RMprotclass)
RMprotclass.columns = ['orgid', 'accession','length','type','cluster','rmT1','rmT2','rmT3','rmT?']
RMprotclass.to_csv('data/RMsearch/RMprotclass.csv',index=False)
