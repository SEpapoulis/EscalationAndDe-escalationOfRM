
import subprocess,tempfile,shutil
import pandas as pd
import numpy as np
from .datastructures import Fasta, ORG

def alignment_conflict(se1,se2,threshold = .75):
    #if the % of the smaller protein is > threshold,
    #we return true

    #se = (start,end)
    #s = start
    #e = end
    #ds = distance from s1 to s2 (aka overhang)
    #de = distance from e1 to e2 (aka overhang)
    s1,e1=se1
    len1=e1-s1
    s2,e2=se2
    len2=e2-s2
    ds=abs(s1-s2)
    de=abs(e1-e2)

    if s1 < s2:
        bs = s1
    else:
        bs = s2
    if e1 > e2:
        be = e1
    else:
        be = e2
    #be = border of end alignment (end coordiante that is the biggest)
    #bs = border of the start alignemnt (start that is the smallest)
    #len_t = total length of area composing both alignments
    #dse = length from total area where ds and de do not overlap (alignment length)
    len_t = be - bs
    dse = len_t - ds - de


    if dse <= 0:
        return(False)
    else:
        if len1 < len2:#find the smaller alignment
            smaller = len1
        else:
            smaller = len2
        #if 75% of the smaller alignment is covered by the bigger alignment, compete them!
        if dse / smaller > threshold:
            return(True)
        else:
            return(False)

def compete_domains(HMMresultsdf,percent_aligned = .75):
    #NOTE: assumes domtblout output
    #removes domains that are out competed by other domains
    #when aligned to the same protein

    key_dat={} #indexing matrix based on query lable
    for line in HMMresultsdf.itertuples(index=False):
        hit=line[0]
        if hit in key_dat:
            key_dat[hit].append(list(line))
        else:
            key_dat[hit]=[list(line)]

    #We will pick which hits to keep by evaluating per protein
    for el in key_dat:
        #There is more than one hit to this sequences, initate competition
        if len(key_dat[el]) > 1:
            out=[]
            eliminated=[]
            for i in range(0,len(key_dat[el])):#for each alignment to el
                align_start=int(key_dat[el][i][17])
                align_end=int(key_dat[el][i][18])
                eval1=float(key_dat[el][i][6])
                a1=(align_start,align_end)
                for i2 in range(i+1,len(key_dat[el])):#compare following hits to look for things to eliminate
                    align2_start=int(key_dat[el][i2][17])
                    align2_end=int(key_dat[el][i2][18])
                    eval2=float(key_dat[el][i2][6])
                    a2=(align2_start,align2_end)
                    if alignment_conflict(a1,a2,threshold = percent_aligned) and i2 not in eliminated:#if we are not comparing the same and there is a conflict
                    #evalue implicitly uses sequence length, threfore:
                    #if sequences overlap in _alignment_conflict, pick the better one
                        if eval1 < eval2:
                            eliminated.append(i2)
                        else:
                            eliminated.append(i)
            for i in range(0,len(key_dat[el])):
                if i not in eliminated:
                    out.append(key_dat[el][i])

            key_dat[el]=out

    #rebuild datastructure
    temp=[]
    for el in key_dat:
        for line in key_dat[el]:
            temp.append(line)
    HMMresultsdf_competed = pd.DataFrame(temp)
    HMMresultsdf_competed.columns = list(HMMresultsdf.columns)
    return(HMMresultsdf_competed)

def collate(BLASTexceptions,probefasta,BLASTresults,HMMresults,fasta,orgid,MINprcnt_aln=.75,FalsePositiveHMM_set=set(),tenative_positiveHMM=set()):

    #BLASTexception fasta holds proteins exceptable to keep with blast evidence only
    BLASTexception=Fasta(BLASTexceptions)
    exception_set = set(BLASTexception.get_accessionlist())
    
    
    #delets alignments not meeting % aligned criteria, then mapping subjects with best hit
    if not BLASTresults.empty:
        subject_besthit={}
        BLASTresults['query_len']=[probefasta.get_length(prot) for prot in BLASTresults['query']]
        blastM = BLASTresults[BLASTresults["alignment length"]/BLASTresults['query_len'] >= MINprcnt_aln]
        #blast files are ordered with the best hit for each query first. iterate and assign best hit to subject
        for row in blastM.itertuples(index=False):
            query,subject,evalue,bitscore=row[0],row[1],row[-3],row[-2]
            if subject in subject_besthit:
                #if current row has lower evalue, assign
                if subject_besthit[subject][1] > evalue:
                    subject_besthit[subject]= (query,evalue,bitscore)
                #if evalues are equal, defer to bitscore (not sure this happens as evalue is calcualted with bitscore?)
                elif subject_besthit[subject][1] == evalue and subject_besthit[subject][2] < bitscore:
                    subject_besthit[subject]= (query,evalue,bitscore)

            else:
                subject_besthit[subject]= (query,evalue,bitscore)

    else:
        blastM=BLASTresults
        subject_besthit={}




    #hmmsearch results using hmmsearchspace, domains competed and poor alignments removed
    if not HMMresults.empty:
        HMMresults = compete_domains(HMMresults,MINprcnt_aln)#compete


    #building seqeunce length and product dictionaries for later
    accession_seqlength={}
    accession_product={}
    for header in fasta:
        accession=header.split(' ')[0][1:]
        accession_seqlength[accession]=len(fasta[header])
        accession_product[accession]='N/A'
        for el in header.split(' ')[1:]:
            if 'product=' in el:
                product=el.strip('[]').split('=')[-1]
                accession_product[accession] = product
        
    ##############################################################

    #adding subjects that were detected with blast and hmmer
    
    if not blastM.empty:
        hitaccessions=set(blastM['subject'])
    else:
        hitaccessions=set()

    #mapping proteins to domains and adding prots to hit set
    accession_domains={}
    for row in HMMresults.itertuples(index=False):
        prot = row[0]
        hitaccessions.add(prot)
        domain = row[3]
        if prot in accession_domains:
            accession_domains[prot].append(domain)
        else:
            accession_domains[prot]=[domain]

    #organzing based on gen bank locus, then sorting based on location
    #fasta files are built in order, thus, use it as a reference
    locus_hitaccessions={}
    for el in hitaccessions:
        locus=el.split('|')[0]
        if locus in locus_hitaccessions:
            locus_hitaccessions[locus].append(el)
        else:
            locus_hitaccessions[locus]=[el]
    ordered_genelist=[]
    for header in fasta:
        accession=header.split(' ')[0][1:]
        if accession in hitaccessions:
            ordered_genelist.append(accession)


    template=[]
    #collect all information for output file
    for accession in ordered_genelist:
        add_seq=False
        accession=accession.split()[0]
        try:
            domains=accession_domains[accession]
        except KeyError:
            domains='N/A'
        length=str(accession_seqlength[accession])
        product=accession_product[accession]
        try:
            topBLASThit_dat=subject_besthit[accession]
            topBLASThit=';'.join([str(el) for el in topBLASThit_dat])
        except KeyError:
            topBLASThit='N/A'

        #if there are domains, add to final output if all domains are not false positives
        if domains!='N/A':
            falsepos=False
            tenative=False
            positive=False
            #if there is a domain of interest, add it
            #if there are false positves with domains of interest, keep it
            #if there are
            for domain in domains:
                #domain is not a false positive or all domains
                if domain in FalsePositiveHMM_set:
                    falsepos=True
                elif domain in tenative_positiveHMM:
                    tenative=True
                else:
                    positive=True
            domains = ';'.join(domains)
            if positive:
                add_seq=True
            elif tenative and not falsepos:
                add_seq=True
        #if there are not domains detected, but was hit with exception list, add to final output
        if topBLASThit.split(';')[0] in exception_set:
            add_seq=True
        if add_seq:
            template.append([accession,domains,length,product,topBLASThit])
    return(template)


def format_hmmsearchfile(raw_inputfile):
    #turns hmmsearch output into tabdelimited file
    raw_hmmresults=[]
    with open(raw_inputfile,'r') as f:
        for line in f:
            if line[0] != '#':
                temp=line.strip()
                temp=temp.split(' ')
                raw_hmmresults.append(temp)
    hmmresults=[]
    for line in raw_hmmresults:
        temp=[]
        for el in line:
            if el != '':
                temp.append(el)
        newline = temp[0:22]
        description=''
        for el in temp[22:]:
            description=description+' '+el
        newline=newline+[description[1:]]
        hmmresults.append(newline)
    return(hmmresults)

def Comprehensive_Scan(target_dir, results_dir, probelibrary, HMM,
BLASTexceptions, threads=1, evalue='.00001',MINprcnt_aln=.75, MAXprcnt_aln=1.25,FalsePositiveHMMs='',TenativePositiveHMMs_list=[]):
    #Comprehensive_Scan assumes that the first element in file.split('.') is the id of the organism
    #Type Conversion
    evalue=str(evalue)
    threads=str(threads)

    #checking files that have already been completed to not compute against
    #note: dir lists drectory in both linux and windows
    possibletarget_files=str(subprocess.check_output(['dir', target_dir]),'utf-8').split()

    result_files=str(subprocess.check_output(['dir', results_dir]),'utf-8').split()

    #Checking if some organisms have been already evaluated
    idcomplete=[]
    for file in result_files:
        if 'bioscan.csv' in file:
            orgid=file.split('.')[0]
            idcomplete.append(orgid)
    try:
        with open(results_dir+'No_element_orgs.txt',"r") as f:
            for line in f:
                idcomplete.append(line.strip())
    except FileNotFoundError:
        pass
    target_files=[]
    for target_file in possibletarget_files:
        orgid=target_file.split('.')[0]
        if orgid not in idcomplete and ('.gbff' in target_file or '.fasta' in target_file):
            target_files.append(target_file)

    #Gathering false positive HMMS
    false_posHMMs=set()
    if FalsePositiveHMMs:
        f = open(FalsePositiveHMMs,'r')
        for line in f:
            false_posHMMs.add(line.strip())
        f.close()

    #defining columns for blast and HMMer

    hmmer_col=[ 'target name','target accession','tlen','query name','accession','qlen',
    'E-value','score','bias','#','of','c-Evalue','i-Evalue','score','bias','from (hmm coord)',
    'to (hmm coord)','from (ali coord)','to (ali coord)','from (env coord)','to (env coord)',
    'acc','description of target'
    ]
    blast_col=["query","subject","percent identity","alignment length",
            "mismatches","gap openings","query start","query end","subject start",
            "subject end","evalue","bitscore"]
    results_col=["accession","Domain(s)Found","length","product","TopBLAST;evalue;bitscore"]

    #book keeping lists
    failed_to_launch=open(results_dir+'failed.txt',"a+") #list of ids
    genome_stats=open(results_dir+'orgid_dat.csv',"a+")#orgid,bp,num_contigs
    element_fasta=open(results_dir+'ElementsFound.fasta',"a+")
    noRM_orgs=open(results_dir+'No_element_orgs.txt',"a+")

    #Launching pipeline for each target file
    i_total = len(target_files)
    for i,target_file in enumerate(target_files):
        orgid=target_file.split('.')[0]
        print('Analyzing: '+orgid)

        #1) Create a temp directory for IO operations
        tempdir=tempfile.mkdtemp(prefix=orgid+'_',suffix='_BioScan_temp')
        
        #2) write fasta
        if '.gbff' in target_file:
            org=ORG(target_dir+target_file)
            fastaCDS=org.get_CDSfasta(include_pseduo=True,gencode=11)
        elif '.fasta' in target_file:
            fastaCDS=Fasta(target_dir+target_file)
        
        if not fastaCDS:
            orgid=target_file.split('.')[0]
            failed_to_launch.write(orgid+'\n')
            print('Cannot find coding sequences, Exiting')

        #Fasta file created, continue with pipeline            
        else:
            fastafile=tempdir+'/'+orgid+'.fasta'
            fastaCDS.write_file(fastafile)

            #3) make blast database to search with threading
            dbfile=tempdir+'/BLASTdb_'+orgid
            db_arg=['makeblastdb','-in',fastafile,'-dbtype','prot','-out',dbfile]
            subprocess.call(db_arg)

            #3)Blast probelibrary against database
            blast_results=tempdir+'/'+orgid+'.blastp.csv'
            blastp_arg=['blastp','-query',probelibrary, '-db',dbfile,'-evalue', evalue,
                         '-num_threads', threads, '-outfmt', '6', '-out', blast_results]
            subprocess.call(blastp_arg)
            try:
                blastresults = pd.read_csv(blast_results,header=None,sep='\t')
                blastresults.columns=blast_col
            except pd.errors.EmptyDataError:
                blastresults=pd.DataFrame(columns=blast_col)
                

            #4)Hmmsearch
            hmmsearch_resultsfile=tempdir+'/'+orgid+'.hmmout'
            hmmlog=tempdir+orgid+'hmmsearch.log'
            hmm_arg=['hmmsearch','--cut_ga', '--domtblout', hmmsearch_resultsfile, HMM, fastafile]
            log=subprocess.check_output(hmm_arg)
            hmmresults = pd.DataFrame(format_hmmsearchfile(hmmsearch_resultsfile))
            if not hmmresults.empty:
                hmmresults.columns=hmmer_col


            #5)collate blast and hmmer results
            probe_fasta= Fasta(probelibrary)
            search_results=collate(probefasta=probe_fasta,BLASTexceptions=BLASTexceptions, BLASTresults=blastresults, HMMresults=hmmresults,fasta=fastaCDS,
                                    orgid=orgid,MINprcnt_aln=MINprcnt_aln,FalsePositiveHMM_set=false_posHMMs, tenative_positiveHMM=set(TenativePositiveHMMs_list))
            search_results = pd.DataFrame(search_results)

            #5) Saving results and organims metadata

            #search results
            if not search_results.empty:
                search_results.columns = results_col
                search_results.to_csv(results_dir+orgid+'.bioscan.csv',index=False)
            else:
                noRM_orgs.write(orgid+'\n')
            
            #saving metadata
            DNAfasta=org.get_DNAfasta()
            genomesize=0
            for dna in DNAfasta:
                genomesize+=len(DNAfasta[dna])
            genome_stats.write(','.join([orgid,str(genomesize),str(len(DNAfasta))])+'\n')

            #saving proteins for future use
            for row in search_results.itertuples(index=False):
                accession,doms,length,product,topblast=list(row)
                if doms == 'N/A':
                    evidence = topblast
                else:
                    evidence = doms
                seq=fastaCDS[accession]
                element_fasta.write(' '.join(['>'+accession, 'alignment_evidence='+evidence, 'orgid='+orgid,
                'length='+length,'annotated_product='+product])+'\n')
                element_fasta.write(seq+'\n')
            print(orgid+' Complete {}/{}\n########################\n\n'.format(str(i+1),str(i_total)))

            #6) clean temp directory to avoid filling up hdd
            shutil.rmtree(tempdir)
    #closing operations
    noRM_orgs.close()
    element_fasta.close()
    genome_stats.close()
    failed_to_launch.close()