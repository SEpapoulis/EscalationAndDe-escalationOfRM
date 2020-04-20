# Resource availability and viral DNA methylation drive the diversity and abundance of Restriction Modification Systems

Bioinformatic and Mathematical Modeling Results

### Dependencies
Python 3 (>=3.7)  
BLAST >=2.9.0+  
HMMER >= 3.1b2  
[cdhit](https://github.com/weizhongli/cdhit)  
[hhsuite](https://github.com/soedinglab/hh-suite)  

Bioinformatic analysis is limited to OS/X and Linux due to alignment software limitations, while numerical simulations are OS independent as they are completed through [SciPy](https://docs.scipy.org/doc/scipy/reference/integrate.html)


## Installation

Navigate to your directory of choice and clone our github repo  
`https://github.com/SEpapoulis/EscalationAndDe-escalationOfRM.git`  
From here, launch any of our jupyter notebooks to view our code, or look at demo/Demo.ipynb for a easily executable sample of our analysis!

#### Easy Installation of Dependencies
We __highly__ recommend installing the [Anaconda data science package](https://www.anaconda.com/distribution/) as all software/modules required is availble through conda install. After anaconda is installed, enter the following in the conda prompt:
`conda install -c bioconda hmmer`  
`conda install -c bioconda blast`  
`conda install -c bioconda cd-hit`  
`conda install -c bioconda hhsuite`  

#### Expert Installation of Dependencies
Independently install all of the dependencies on your machine.
[cdhit](https://github.com/weizhongli/cdhit)  
[hhsuite](https://github.com/soedinglab/hh-suite)  
[HMMER](http://hmmer.org/download.html)  
[BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  

## Data Demo
We have provided a demo of RM pipeline as well as some sample code to run numerical simulations using our memory model. Open the jupyer notebook and start executing cells! The notebook highlights modifications you can make to run other assemblies from Refseq. While the sample bioinformatic analysis with *Microcystis* takes ~15min, numerical simulations may take between one and two hours depending on your machine. If you wish to see all notebooks used in this research, please see the NotebooksAndData directory.


### Viewing project code
This code documents the results used to deepen our understanding of the selective pressures that govern the loss and gain of Restriction Modification Systems among prokaryotes. This project has been documented in exacutable [Jupyter Notebooks](https://jupyter.org/) and all imports are listed at the begining of notebook, with the exception of R code where libraries are imported per script. Jupyter Notebooks can be rendered in github, however, if notebooks fail to render, they can be downloaded and viewed locally. To view Jupyter Notebooks locally, simply install the [Anaconda data science package](https://www.anaconda.com/distribution/). 


### Project Files
Code used to generate main body figures can be viewed in the Manuscript_Figures.ipynb while supplementary figures can be viewed in the Manuscript_Figures_SI.ipynb notebook, except for figures S7-S11, which can be found in NotebooksAndData/RM_Database/definitions/Defining_RM_types_5-3-2018.ipynb. Bellow, we describe the contents of our project folders.

1. RM_Database - All database files needed to recapitulate RM annotation
   1. BLASTexceptions.fasta - Sequences without HMMs usef to find RM genes via BLAST
   2. Falsepos_HMMs.txt - HMMs that covaried with false positives (mostly helicases)
   3. RM_HMMs.txt - HMMs found in non-puatative Restriction Modification Sysetms
   4. non-putative_rebase.fasta - A reformatted file of [NEB's non-putative protein sequences](ftp://ftp.neb.com/pub/rebase/protein_mini_reg_seqs.txt)
   5. definitions (directory) - Our initial inquery used to define our database. SI figures are found here in Defining_RM_types_5-3-2018.ipynb



2. RMsearch
   1. src (directory) - holds source code for RM searches
   2. RMsearch_12-11-2018.ipynb - Notebook for RM annotation
   3. data (directory) - outputs of RMsearch notebook


3. GenMemODE
   1. src (directory) - holds source code needed for all numerical simulations
   2. Computational_Figures.ipynb - initial figures used for understanding our simulations
   3. GenMemODE.ipynb - Numerical Simulations



#### Additional Database files for hhsuite and HMMs used in this study
[uniprot](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/uniprot20_2016_02.tgz)  
[protein databank](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/pdb70_14Sep16.tgz)  
[pfam 31.0](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz)  


