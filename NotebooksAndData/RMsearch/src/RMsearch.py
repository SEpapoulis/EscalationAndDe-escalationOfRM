import sys

sys.path.append('/home/zinserlab/bin/Bio_workbench')

from Bioscan import Search_Strategies
from datetime import datetime

##path variables##
db_location='/media/zinserlab/4TbWesterDigital/Databases/RM/'
gb_dir='/media/zinserlab/4TbWesterDigital/Prokaryotes/'
probelibrary= db_location+'non-putative_rebase.fasta'
HMM_searchspace= db_location+'RMwithFalsepos.hmm'
BLASTexceptions= db_location+'BLASTexceptions.fasta'
falsepos=db_location+'Falsepos_HMMs.txt'
results= '/home/zinserlab/Documents/Spiro/projects/RMsearch/data/search_results/'
startTime = datetime.now()

#when probe library is blast exceptions, 10x faster
Search_Strategies.Comprehensive_Scan(target_dir=gb_dir,results_dir=results,probelibrary=BLASTexceptions,
               HMM=HMM_searchspace,BLASTexceptions=BLASTexceptions,threads=3,FalsePositiveHMMs=falsepos,TenativePositiveHMMs_list=['ResIII'])


print ('runtime: '+str(datetime.now() - startTime))
