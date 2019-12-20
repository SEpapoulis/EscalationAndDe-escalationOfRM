from subprocess import check_output
import gzip
Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
Base2 = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
Base3 = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

gencodedat=[
            ['11',
            'Bacterial, Archaeal and Plant Plastid',
            'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
            '---M------**--*----M------------MMMM---------------M------------']]

amb_baseset={'R':set(['A','G']),
            'Y':set(['C','T']),
            'S':set(['G','C']),
            'W':set(['A','T']),
            'K':set(['G','T']),
            'M':set(['A','C']),
            'B':set(['C','T','G']),
            'D':set(['A','G','T']),
            'H':set(['A','C','T']),
            'V':set(['A','C','G']),
            'N':set(['A','C','T','G']),
            }

def is_ambiguous(nt):
    if nt.upper() in amb_baseset:
        return(True)
    return(False)

def has_ambiguous(codon):
    if is_ambiguous(codon[0]) or is_ambiguous(codon[1]) or is_ambiguous(codon[2]):
        return(True)
    return(False)

#returns alternate codon forms
#if there isnt an alternate codon,
#returns empty list
def get_alternate_codons(codon):
    alts=[]
    p1=codon[0]
    p2=codon[1]
    p3=codon[2]
    if is_ambiguous(p1):
        p1 = amb_baseset[p1]
    if is_ambiguous(p2):
        p2 = amb_baseset[p2]
    if is_ambiguous(p3):
        p3 = amb_baseset[p3]
    for el1 in p1:
        for el2 in p2:
            for el3 in p3:
                alt = el1+el2+el3
                if alt != codon:
                    alts.append(alt)
    return(alts)

def load_transtbls():
    transtbl_codondicts={}
    
    for dat in gencodedat:
        tblid=int(dat[0])
        AAs=dat[2]
        ss=dat[3]
        codon_AA={}
        codon_ss={} #start/stop
        
        for i,AA in enumerate(AAs):
            codon = ''.join([Base1[i],Base2[i],Base3[i]])
            codon_AA[codon] = AA
            codon_ss[codon] = ss[i]
        transtbl_codondicts[tblid] = (codon_AA,codon_ss)
    return (transtbl_codondicts)
    

def load_NTcomplements():
    nt_complement={
        'a':'t','t':'a','g':'c','c':'g', 
        'n':'n','r':'y','y':'r','s':'s',
        'w':'w','k':'m','m':'k','b':'v',
        'v':'b','d':'h','h':'d'}
    return(nt_complement)

transtbls = load_transtbls()
comp=load_NTcomplements()


class Fasta():
    def __init__(self,file=''):
        self.headerlist=[] #will retain order of things added to the dictionary
        self.accession_headerindex={} #maps accession to full header
        self.header_seq={} #full header to sequence

        if file:
            self.load_file(file)

    def __setitem__(self, header, seq):
        #if user does not add > in header, do it for them
        if header[0] != '>':
            header = '>'+header
        self.headerlist.append(header)
        self.header_seq[header]=seq
        accession = header[1:].split()[0] #removing '>' then isolating accession
        self.accession_headerindex[accession]=len(self.headerlist)-1
        
    
    def __getitem__(self, key):
        try:
            seq=self.header_seq[key] #are they getting the sequence with the full header?
        except KeyError: #they must be searching by accession
            headerindex=self.accession_headerindex[key]
            seq=self.header_seq[self.headerlist[headerindex]]
        return(seq)

    def __iter__(self):
        return(iter(self.headerlist))

    def __len__(self):
        return(len(self.headerlist))


    def get_header(self,accession):
        i = self.accession_headerindex[accession]
        return(self.headerlist[i])

    def get_length(self,accession):
        header = self.get_header(accession)
        return(len(self[header]))

    def get_sequence(self,accession):
        header = self.get_header(accession)
        return(self[header])

    def get_accessionlist(self,retrieve=set()):
        temp=[]
        for el in self.headerlist:
            accession=el.split()[0]
            if retrieve and accession in retrieve:
                temp.append(accession[1:])
            else:
                temp.append(accession[1:])
        return(temp)

    def load_file(self,fastafile):
        f=open(fastafile,'r')
        seq=""
        key=''
        for line in f:
            if line[0] == '>':
                if key!='':
                    self[key]=seq
                key=line.strip()
                seq=""
            else:
                temp=line.strip()#removing end of line
                temp=temp.replace(" ","")#removing any possible spaces
                seq+=temp
        if key:#just incase file is empty
            self[key]=seq#write last seq
        f.close()

    def write_file(self,filename):#write fasta stored in dictionary to file
        f=open(filename,'w')
        for key in self:
            f.write(key+'\n')
            dat=self[key]
            for x in range(0,len(dat),100):
                if x+100 <= len(dat):
                    f.write(dat[x:x+100]+'\n')
                else:
                    f.write(dat[x:len(dat)]+'\n')
        f.close()


class ORG:
    def __init__(self,file,taxid=''):
        self.file=file
        self.GENEtable=[]       #[0] locus_location, rest are data feilds
        self.CDStable=[]        #[0] locus_location, rest are data feilds
        self.rRNAtable=[]       #[0] locus_location, rest are data feilds
        self.tRNAtable=[]       #[0] locus_location, rest are data feilds
        self.SOURCEtable=[]     #[0] locus_location, rest are data feilds
        self.misc_MISCtable={}  #stores unspecified tables in dictionary
        self.locus_DNA={}
        self.organism=''
        self.strain=''
        if taxid:
            self.taxid=taxid
        else:
            self.taxid=''

        if self.file.split('.')[-1] == 'gbff' or (self.file.split('.')[-1] == 'gz' and self.file.split('.')[-2] == 'gbff' ):
            try:
                self._ParseGenbank()
                for line in self.SOURCEtable:
                    for dat in line:
                        if 'organism' in dat:
                            self.organism=dat.split('=')[1]
                        elif 'db_xref=taxon:' in dat and not taxid:
                            self.taxid=dat.split(':')[1]
                        elif 'strain=' in dat:
                            self.strain=dat.split('=')[1]
                    if self.organism and self.strain and self.taxid:
                        break
            except EOFError:
                print("EOF not found in "+file+",please redownload")

    def _process_FEATURES(self,FEATURESstr,current_locus):
        dat=FEATURESstr.split('\t')
        dat.pop(0)#removing empty line
        for el in dat:
            fields=el.replace('"','').split('/')
            field_type,location=fields.pop(0).split()
            location = current_locus + '|' + location
            fields.insert(0,location)
            if field_type == 'gene':
                self.GENEtable.append(fields)
            elif field_type == 'CDS':
                self.CDStable.append(fields)
            elif field_type == 'rRNA':
                self.rRNAtable.append(fields)
            elif field_type == 'tRNA':
                self.tRNAtable.append(fields)
            elif field_type == 'source':
                self.SOURCEtable.append(fields)
            else:
                try:
                    self.misc_MISCtable[field_type].append(fields)
                except KeyError:
                    self.misc_MISCtable[field_type]=[fields]

    def _ParseGenbank(self):
        FEATURES=[] # cds,gene,rRNA,ect...
        _dnatemp=[] #holds dna of locus being reading
        current_locus=''
        read_features=False
        read_origin=False
        gz=False
        if self.file.split('.')[-1] == 'gz':
            gz=True
            gbff=gzip.open(self.file)
        else:
            gbff=open(self.file,'r')

        #readfile#
        for line in gbff:
            if gz:
                line=str(line,'utf-8')

            #locus found, write data if not first locus and reset
            if 'LOCUS' == line[:5]:
                if current_locus:
                    self.locus_DNA[current_locus]=''.join(_dnatemp)
                    self.temp=''.join(FEATURES)
                    self._process_FEATURES(''.join(FEATURES),current_locus)
                current_locus=line.split()[1]
                FEATURES=[]
                _dnatemp=[]
                read_features=False
                read_origin=False

            #dna is next#
            elif 'ORIGIN' == line[:6]:
                read_origin = True
                read_features = False

            #Features are next#
            elif 'FEATURES' == line[:8]:
                read_features = True

            #now reading features, avoiding CONTIG line#
            elif read_features:
                #member of current feature?
                if '                     ' == line[:21]:
                    line = line.strip()
                    FEATURES.append(line)
                #new feature, add \t for split later
                elif '     ' == line[:5]:
                    line = line.strip()
                    FEATURES.append('\t'+line)
                else:
                    read_features = False

            #now reading DNA#
            elif read_origin:
                if '//' == line[:2]:
                    read_origin = False
                else:
                    line=line.strip()
                    for el in line.split(' ')[1:]:
                        _dnatemp.append(el)
        #file iteration has ended, write last data
        self.locus_DNA[current_locus]=''.join(_dnatemp)
        self._process_FEATURES(''.join(FEATURES),current_locus)
        gbff.close()

    def reverse_complement(self,DNA):
        complement = [comp[nt] for nt in DNA.lower()]
        return(''.join(complement[::-1]).upper())


#returns DNA 5' to 3'
#this should get simplified?
    def get_DNA(self,locus_location,locus=''):
        if locus:
            location=locus_location
        else:
            locus,location=locus_location.split('|')[0],locus_location.split('|')[1]
            location = location.replace(">", "") #removing truncation marker
            location = location.replace("<", "") #removing truncation marker
        if location[0:4] == 'join':
            location=location[5:-1]
            loclist=location.split(',')
            DNA=''
            for el in loclist:
                DNA=DNA+self.get_DNA(el,locus)
            return(DNA)
        elif location[0:4] == 'comp':
            location=location[11:-1]
            DNA=self.get_DNA(location,locus)
            rcDNA=self.reverse_complement(DNA)
            return(rcDNA)
        else:
            if '..' not in location:#condition for single nt!
                DNA=self.locus_DNA[locus][int(location)-1]
                return(DNA)
            l,r=location.split('..')
            DNA=self.locus_DNA[locus][int(l)-1:int(r)]   #genbank starts array at 1... subtracting 1 to align nt with python
            return(DNA)

#requires revision
    def get_DNAfasta(self):
        fasta=Fasta()
        for locus in self.locus_DNA:
            header='>'+locus
            fasta[header]=self.locus_DNA[locus]
        return(fasta)
            

    def get_pseudogenes(self):
        fasta=Fasta()
        for line in self.GENEtable:
            for el in line[1:]:
                if 'pseudo' in el.lower():
                    fasta['>'+line[0]] = self.get_DNA(line[0])
                    break
        return(fasta)

    def translate(self,locus, DNA, gencode, start_index=0, stop_index=0, iscomplement=False, DNAcircular = False):
        DNA=DNA.upper()
        possible_start=[]
        start=start_index
        AA=[]
        codon_AA,codon_SS=transtbls[gencode]
        fasta = Fasta()
        if stop_index:
            stop_index = stop_index - (stop_index%3)
        else:
            stop_index = (len(DNA) - start_index) - ((len(DNA)- start_index)%3)
        for i in range(start_index,stop_index,3):
            codon=DNA[i:i+3]
            if has_ambiguous(codon):
                alternates = get_alternate_codons(codon)
                possible_AA = [codon_AA[c] for c in alternates]
                if set(possible_AA) == 1:
                    AA.append(possible_AA[0])
                else:
                    AA.append('X')

                possible_ss = [codon_SS[c] for c in alternates]
                #if degeneracy leads to all codons still being start, add as possible start
                if set(possible_ss) == set(['M']):
                    possible_start.append(str(i+1))
            else:
                AA.append(codon_AA[codon])
                if codon_SS[codon] == 'M':
                    possible_start.append(str(i+1))

            #Check for stop codon to add to fasta
            if AA[-1] == '*':
                if AA[0] != '*': #checking to make sure the first AA added was not stop
                    header = ['>',locus,'|',str(start+1),'..',str(i),' [possible_starts:'+','.join(possible_start)+']']
                    if iscomplement:
                        header.insert(6,')')
                        header.insert(3,'complement(')
                    fasta[''.join(header)] = ''.join(AA[:-1])
                AA=[]
                possible_start=[]
                start = i+3
        #for the last sequence
        header = ['>',locus,'|',' [possible_starts:'+','.join(possible_start)+']']
        
        #special case where AA is empty
        if not AA:
            location=''
        #we want to include the last AA because it is not a *
        elif AA[-1] == '*':
            seq=AA[:-1]
            location =[str(start+1),'..',str(i)]
        else:
            seq=AA
            location =[str(start+1),'..',str(i+3)]
        header.insert(3,''.join(location))
        if iscomplement:
            header.insert(4,')')
            header.insert(3,'complement(')
        #only add if a location was generated
        if location:
            fasta[''.join(header)] = ''.join(seq)
        return(fasta)


#ignores pseudo genes
#this is a fairly old function, should be updated
    def get_CDSfasta(self,product_keywords=list(),include_pseduo=False,gencode=11):
        CDSfasta=Fasta()
        keep = bool
        for dat in self.CDStable:
            if 'translation=' in dat[-1]:
                product = 'NA'
                protein_id = 'NA'
                for el in reversed(dat): #usually at the end, should increase search speed
                    if 'product=' in el:
                        product=el.replace(" ","_") #removing any spaces from product
                    elif 'protein_id=' in el:
                        protein_id=el
                    elif protein_id != 'NA' and product != 'NA':
                        break
                #>locus|Location [product] [protein_id] [taxid]
                if product_keywords:
                    keep=False
                    for keyword in product_keywords:
                        if keyword.lower() in product.lower(): #making sure it is not case sensitive
                            keep = True
                            break
                else:
                    keep = True
                if keep:
                    header='>'+dat[0]+' ['+product+']'+' ['+protein_id+']'+' [taxid='+self.taxid+']'
                    seq=dat[-1].split('=')[-1]
                    CDSfasta[header]=seq
            elif include_pseduo:
                locus = dat[0].split('|')[0]
                DNA=self.get_DNA(dat[0])
                fasta = self.translate(locus,DNA,gencode)
                seq=[]
                product = 'NA'
                for header in fasta:
                    seq.append(fasta[header])
                for el in reversed(dat): #usually at the end, should increase search speed
                    if 'product=' in el:
                        product=el.replace(" ","_") #removing any spaces from product
                    elif product != 'NA':
                        break
                CDSfasta['>'+dat[0]+' ['+product+']'] = ''.join(seq)
                
        return(CDSfasta)

#requires revision
    #will return a dictionary containing possible 16S of organism
    def get_16Sfasta(self):
        names=['locus','organism=','strain=','taxid=']
        Fasta16S=Fasta()
        #collect all 16S rRNA
        for line in self.rRNAtable:
            for dat in line:
                if '16S' in dat and 'product=' in dat:
                    info=[line[0],self.organism,self.strain,self.taxid]
                    DNA=self.get_DNA(info[0])
                    header='>'+line[0]+' ['+dat+']'
                    for x in range(1,len(info)):
                        el=info[x]
                        if el:
                            header=header+' ['+names[x]+str(el)+']'
                    Fasta16S[header]=DNA
        return(Fasta16S)

    def get_taxid(self):
        return(self.taxid)
    def get_locusDNA(self):
        return(self.locus_DNA)

    def get_locuslist(self):
        temp=[]
        for el in self.locus_DNA:
            temp.append(el)
        return(temp)

    def get_CDStable(self):
        return(self.CDStable)
    def get_GENEtable(self):
        return(self.GENEtable)
    def get_rRNAtable(self):
        return(self.rRNAtable)
